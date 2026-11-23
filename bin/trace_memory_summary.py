#!/usr/bin/env python3
"""
Rank a TCRtoolkit run's memory use by process and by patient.

Joins the Nextflow trace (pipeline_info/execution_trace_*.txt, written on every run) to
pipeline_info/tcrdist_metrics.tsv (one row per patient from TCRDIST3_MATRIX: clone counts,
sparse nnz, radius, wall time, peak RSS) and, when the work directory is local, to the
OBJECT_LOAD lines the single-cell steps log when they load a Seurat / H5AD object.

    trace_memory_summary.py --trace out/pipeline_info/execution_trace_*.txt \
        --tcrdist-metrics out/pipeline_info/tcrdist_metrics.tsv \
        --scan-logs --outdir out/pipeline_info/memory_summary

Writes:
    memory_by_process.tsv  one row per process: tasks, failures, peak RSS (max / median),
                           requested memory, longest realtime, the task with the max
    memory_by_patient.tsv  one row per patient x per-patient process (task tag = patient),
                           with the TCRdist3 clone counts, so memory can be read against
                           n_patient_clones; sorted by peak RSS
    tcrdist_scaling.tsv    TCRDIST3_MATRIX only: peak RSS per clone and per stored pair,
                           to extrapolate the memory a larger patient will need
    objects_loaded.tsv     the Seurat / H5AD objects loaded per task (with --scan-logs)

Several trace files (resumed runs) can be given: the last attempt of each task wins.
"""

import argparse
import glob
import os
import re
import sys

import pandas as pd

UNITS = {"B": 1, "KB": 1024, "MB": 1024 ** 2, "GB": 1024 ** 3, "TB": 1024 ** 4}


def to_bytes(v):
    """'1.2 GB' / '512 MB' / raw byte counts -> bytes; '-' or empty -> NaN."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return float("nan")
    s = str(v).strip()
    if s in ("", "-"):
        return float("nan")
    if re.fullmatch(r"\d+(\.\d+)?", s):
        return float(s)
    m = re.fullmatch(r"([\d.]+)\s*([KMGT]?B)", s, re.I)
    return float(m.group(1)) * UNITS[m.group(2).upper()] if m else float("nan")


def to_seconds(v):
    """'1h 2m 3s' / '4.5s' / '850ms' / raw milliseconds -> seconds."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return float("nan")
    s = str(v).strip()
    if s in ("", "-"):
        return float("nan")
    if re.fullmatch(r"\d+", s):
        return float(s) / 1000.0
    total, found = 0.0, False
    for num, unit in re.findall(r"([\d.]+)\s*(ms|d|h|m|s)", s):
        found = True
        total += float(num) * {"ms": 0.001, "s": 1, "m": 60, "h": 3600, "d": 86400}[unit]
    return total if found else float("nan")


def gb(x):
    return round(x / 1024 ** 3, 2) if pd.notna(x) else float("nan")


def read_traces(paths):
    frames = []
    for i, p in enumerate(paths):
        t = pd.read_csv(p, sep="\t", dtype=str)
        t["_file_order"] = i
        frames.append(t)
    t = pd.concat(frames, ignore_index=True)
    # Traces written with the default fields have no process / tag column: both are in the
    # task name, "WORKFLOW:...:PROCESS (tag)".
    name = t["name"].astype(str)
    if "process" not in t:
        t["process"] = name.str.replace(r"\s*\(.*\)$", "", regex=True)
    if "tag" not in t:
        t["tag"] = name.str.extract(r"\((.*)\)$", expand=False)
    t["process_short"] = t["process"].astype(str).str.split(":").str[-1]
    nan = pd.Series(float("nan"), index=t.index)
    t["peak_rss_b"] = t["peak_rss"].map(to_bytes) if "peak_rss" in t else nan
    t["memory_b"] = t["memory"].map(to_bytes) if "memory" in t else nan
    t["realtime_s"] = t["realtime"].map(to_seconds) if "realtime" in t else nan
    t["attempt_n"] = pd.to_numeric(t["attempt"], errors="coerce").fillna(1) if "attempt" in t else 1
    # Last attempt of each task (a task resumed across trace files keeps its hash).
    key = "hash" if "hash" in t else "name"
    t = t.sort_values(["_file_order", "attempt_n"]).drop_duplicates(key, keep="last")
    return t


def by_process(t):
    rows = []
    for proc, g in t.groupby("process_short"):
        top = g.loc[g["peak_rss_b"].idxmax()] if g["peak_rss_b"].notna().any() else g.iloc[0]
        rows.append({
            "process": proc,
            "tasks": len(g),
            "failed": int(g["status"].isin(["FAILED", "ABORTED"]).sum()) if "status" in g else 0,
            "max_peak_rss_gb": gb(g["peak_rss_b"].max()),
            "median_peak_rss_gb": gb(g["peak_rss_b"].dropna().median()) if g["peak_rss_b"].notna().any() else float("nan"),
            "max_requested_gb": gb(g["memory_b"].max()),
            "max_realtime_h": round(g["realtime_s"].max() / 3600, 2) if g["realtime_s"].notna().any() else float("nan"),
            "max_task": top.get("name", ""),
            "max_task_workdir": top.get("workdir", ""),
        })
    return pd.DataFrame(rows).sort_values("max_peak_rss_gb", ascending=False)


def scan_objects(t):
    """OBJECT_LOAD lines from each task's .command.log (local work directories only)."""
    rows = []
    for _, r in t.iterrows():
        wd = str(r.get("workdir", ""))
        log = os.path.join(wd, ".command.log")
        if not wd or wd.startswith(("s3://", "gs://", "az://")) or not os.path.exists(log):
            continue
        with open(log, errors="replace") as fh:
            for line in fh:
                if "OBJECT_LOAD" not in line:
                    continue
                fields = dict(kv.split("=", 1) for kv in line.strip().split("\t")[1:] if "=" in kv)
                fields.update(task=r.get("name", ""), task_peak_rss_gb=gb(r["peak_rss_b"]))
                rows.append(fields)
    if not rows:
        return None
    d = pd.DataFrame(rows)
    for c in ("file_mb", "object_mb", "cells", "features", "r_max_used_mb", "peak_rss_mb"):
        if c in d:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    return d.sort_values("file_mb", ascending=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--trace", nargs="+", required=True, help="Nextflow trace file(s); globs allowed")
    ap.add_argument("--tcrdist-metrics", default=None, help="pipeline_info/tcrdist_metrics.tsv")
    ap.add_argument("--scan-logs", action="store_true",
                    help="read OBJECT_LOAD lines from each task's local .command.log")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()

    paths = sorted({p for pat in args.trace for p in glob.glob(pat)})
    if not paths:
        sys.exit(f"No trace file matches {args.trace}")
    os.makedirs(args.outdir, exist_ok=True)
    t = read_traces(paths)
    print(f"[trace_memory_summary] {len(t)} tasks from {len(paths)} trace file(s)")

    bp = by_process(t)
    bp.to_csv(os.path.join(args.outdir, "memory_by_process.tsv"), sep="\t", index=False)
    print(bp.head(10).to_string(index=False))

    # Per-patient processes carry the patient as their tag (TCRDIST3_MATRIX, GIANA,
    # GLIPH2, PATIENT_CONCATENATE ...). Other tags (sample, project) are left out.
    metrics = None
    if args.tcrdist_metrics and os.path.exists(args.tcrdist_metrics):
        metrics = pd.read_csv(args.tcrdist_metrics, sep="\t", dtype={"patient": str})
    if metrics is not None and "tag" in t:
        pt = t[t["tag"].isin(set(metrics["patient"]))].copy()
        pt = pt.merge(metrics.rename(columns={"patient": "tag",
                                              "peak_rss_mb": "tcrdist_script_peak_rss_mb",
                                              "wall_time_s": "tcrdist_script_wall_s"}),
                      on="tag", how="left")
        out = pt.assign(patient=pt["tag"], peak_rss_gb=pt["peak_rss_b"].map(gb),
                        requested_gb=pt["memory_b"].map(gb),
                        realtime_h=(pt["realtime_s"] / 3600).round(2))
        keep = ["patient", "process_short", "status", "attempt", "peak_rss_gb", "requested_gb",
                "realtime_h", "n_samples", "n_input_rows", "n_unique_cdr3_beta", "n_patient_clones",
                "radius", "nnz", "n_clusters", "tcrdist_script_peak_rss_mb", "tcrdist_script_wall_s"]
        out = out[[c for c in keep if c in out.columns]].rename(columns={"process_short": "process"})
        out.sort_values("peak_rss_gb", ascending=False).to_csv(
            os.path.join(args.outdir, "memory_by_patient.tsv"), sep="\t", index=False)

        td = out[out["process"] == "TCRDIST3_MATRIX"].copy()
        if len(td):
            td["gb_per_1k_clones"] = (td["peak_rss_gb"] / (td["n_patient_clones"] / 1000)).round(3)
            td["bytes_per_stored_pair"] = (td["peak_rss_gb"] * 1024 ** 3 / td["nnz"].where(td["nnz"] > 0)).round(1)
            td.sort_values("n_patient_clones", ascending=False).to_csv(
                os.path.join(args.outdir, "tcrdist_scaling.tsv"), sep="\t", index=False)
            print("\nTCRDIST3_MATRIX, largest patients:")
            print(td.sort_values("n_patient_clones", ascending=False)
                    [["patient", "n_patient_clones", "nnz", "peak_rss_gb", "realtime_h"]].head(10).to_string(index=False))

    if args.scan_logs:
        objs = scan_objects(t)
        if objs is not None:
            objs.to_csv(os.path.join(args.outdir, "objects_loaded.tsv"), sep="\t", index=False)
            print(f"\n{len(objs)} Seurat / H5AD loads logged; largest:")
            print(objs.head(5).to_string(index=False))
    print(f"[trace_memory_summary] wrote {args.outdir}")


if __name__ == "__main__":
    main()
