#!/usr/bin/env python3
"""
COMPARTMENT_OVERLAP: clonotype sharing between sample compartments within each patient
(e.g. tumor and blood, or any other tissue / body fluid named in the samplesheet).

For every patient sampled in at least two compartments, clone counts are pooled per
compartment (all of the patient's samples in that compartment) and compared pairwise.
The same is done per cohort timepoint (--time-col, e.g. timepoint as binned study days) where a patient was
sampled in two or more compartments at that time.

Inputs : concatenated_cdr3_sorted.tsv (sample, junction_aa, v_call, duplicate_count, ...)
         samplesheet CSV (sample, patient/subject, compartment column, time column)
Outputs: compartment_overlap_summary.tsv  pairwise metrics per patient x stratum
         compartment_flow.tsv             per-compartment frequencies of the clones drawn
                                          in the flow diagram (top-N per compartment) plus
                                          an "Other" remainder row
         compartment_shared_clones.tsv    every clone found in >= 2 compartments of a
                                          patient, with its frequency in each
         compartment_overlap_status.txt   one line: "ok" or why the analysis was skipped
         compartment_samples.tsv          per-sample patient, compartment, timepoint and
                                          cohort time (the Master Summary's sample table)

A run with no compartment column, or no patient sampled in two compartments, still writes
every file (headers only) and states the reason in the status file, so downstream reports
can skip the section cleanly.
"""

import argparse
import itertools
import re
import sys

import numpy as np
import pandas as pd

# Canonical compartment labels; anything unlisted keeps its own (stripped) value.
COMPARTMENT_ALIASES = {
    "blood": "Blood", "pbmc": "Blood", "pbmcs": "Blood", "pb": "Blood",
    "peripheral blood": "Blood", "whole blood": "Blood", "plasma": "Blood",
    "csf": "CSF", "cerebrospinal fluid": "CSF", "cerebrospinal_fluid": "CSF",
    "tumor": "Tumor", "tumour": "Tumor", "til": "Tumor", "tils": "Tumor",
    "tumor tissue": "Tumor", "tumour tissue": "Tumor",
}
# Display / flow order: any other compartment (alphabetical) first, then Blood, then Tumor,
# so blood sits next to every compartment it is compared with.
LAST_COMPARTMENTS = ["Blood", "Tumor"]

PATIENT_CANDIDATES = ("patient", "patient_id", "subject_id", "subject", "donor")

SUMMARY_COLS = [
    "patient", "stratum", "compartment_a", "compartment_b", "samples_a", "samples_b",
    "clones_a", "clones_b", "shared_clones",
    "frac_clones_a_shared", "frac_clones_b_shared",
    "frac_repertoire_a_shared", "frac_repertoire_b_shared",
    "jaccard", "morisita_horn",
]
FLOW_COLS = ["patient", "stratum", "compartment", "clone", "junction_aa", "v_gene",
             "count", "freq", "rank", "shared", "compartments_present"]
SHARED_COLS = ["patient", "stratum", "clone", "junction_aa", "v_gene",
               "compartments_present", "n_compartments"]


def normalise_compartment(x) -> str:
    s = str(x).strip()
    return COMPARTMENT_ALIASES.get(s.lower(), s)


def compartment_sort_key(c: str):
    return (1 + LAST_COMPARTMENTS.index(c) if c in LAST_COMPARTMENTS else 0, c)


def v_gene(v) -> str:
    """Gene-level V call: first call of an ambiguous list, allele dropped (TRBV5-1*01 -> TRBV5-1)."""
    if pd.isna(v):
        return ""
    return str(v).split(",")[0].split("*")[0].strip()


def time_sort_key(t):
    """Numeric part of a time label (d28, 28, Day 28 -> 28); non-numeric labels sort last."""
    m = re.search(r"-?\d+(\.\d+)?", str(t))
    return (0, float(m.group())) if m else (1, str(t))


def morisita_horn(a: np.ndarray, b: np.ndarray) -> float:
    if a.sum() == 0 or b.sum() == 0:
        return float("nan")
    pa, pb = a / a.sum(), b / b.sum()
    den = (pa ** 2).sum() + (pb ** 2).sum()
    return float(2 * (pa * pb).sum() / den) if den > 0 else float("nan")


def write_outputs(summary, flow, shared, status, prefix=""):
    summary.reindex(columns=SUMMARY_COLS).to_csv(f"{prefix}compartment_overlap_summary.tsv", sep="\t", index=False)
    flow.reindex(columns=FLOW_COLS).to_csv(f"{prefix}compartment_flow.tsv", sep="\t", index=False)
    shared.reindex(columns=SHARED_COLS + [c for c in shared.columns if c not in SHARED_COLS]) \
          .to_csv(f"{prefix}compartment_shared_clones.tsv", sep="\t", index=False)
    with open(f"{prefix}compartment_overlap_status.txt", "w") as fh:
        fh.write(status + "\n")
    print(f"[COMPARTMENT_OVERLAP] {status}")


def write_sample_table(meta: pd.DataFrame, args) -> None:
    """Per-sample metadata for the report: written whatever the analysis outcome."""
    def col(name):
        return meta[name] if name and name in meta.columns else pd.Series("", index=meta.index)
    subj = next((c for c in (args.subject_col,) + PATIENT_CANDIDATES if c in meta.columns), None)
    has_time = args.time_col in meta.columns and (meta[args.time_col] != "").any()
    out = pd.DataFrame({
        "sample": meta["sample"],
        "alias": col("alias").where(col("alias") != "", meta["sample"]),
        "patient": col(subj),
        "compartment": col(args.compartment_col).map(lambda x: normalise_compartment(x) if x else ""),
        "timepoint": col(args.fallback_time_col),
        "cohort_time": col(args.time_col) if has_time else col(args.fallback_time_col),
        "cohort_time_col": args.time_col if has_time else args.fallback_time_col,
    })
    out.to_csv("compartment_samples.tsv", sep="\t", index=False)


def analyse_stratum(counts: pd.DataFrame, patient: str, stratum: str, top_n: int,
                    n_samples: dict):
    """counts: clone x compartment count matrix for one patient and stratum."""
    comps = sorted(counts.columns, key=compartment_sort_key)
    counts = counts[comps]
    freqs = counts / counts.sum(axis=0)
    present = counts > 0

    summary = []
    for a, b in itertools.combinations(comps, 2):
        in_a, in_b = present[a], present[b]
        shared = in_a & in_b
        union = (in_a | in_b).sum()
        summary.append({
            "patient": patient, "stratum": stratum, "compartment_a": a, "compartment_b": b,
            "samples_a": n_samples.get(a, 0), "samples_b": n_samples.get(b, 0),
            "clones_a": int(in_a.sum()), "clones_b": int(in_b.sum()),
            "shared_clones": int(shared.sum()),
            "frac_clones_a_shared": shared.sum() / in_a.sum() if in_a.sum() else np.nan,
            "frac_clones_b_shared": shared.sum() / in_b.sum() if in_b.sum() else np.nan,
            "frac_repertoire_a_shared": float(freqs.loc[shared, a].sum()),
            "frac_repertoire_b_shared": float(freqs.loc[shared, b].sum()),
            "jaccard": shared.sum() / union if union else np.nan,
            "morisita_horn": morisita_horn(counts[a].to_numpy(float), counts[b].to_numpy(float)),
        })

    n_present = present.sum(axis=1)
    where = present.apply(lambda r: ";".join(c for c in comps if r[c]), axis=1)

    # Clones drawn in the flow diagram: union of each compartment's top-N.
    drawn = set()
    for c in comps:
        drawn.update(counts[c][counts[c] > 0].sort_values(ascending=False).head(top_n).index)
    flow = []
    for c in comps:
        rank = counts[c].rank(ascending=False, method="first")
        for clone in drawn:
            if counts.at[clone, c] > 0:
                flow.append({"patient": patient, "stratum": stratum, "compartment": c,
                             "clone": clone, "count": int(counts.at[clone, c]),
                             "freq": float(freqs.at[clone, c]), "rank": int(rank[clone]),
                             "shared": bool(n_present[clone] > 1),
                             "compartments_present": where[clone]})
        rest = counts.loc[~counts.index.isin(drawn), c]
        flow.append({"patient": patient, "stratum": stratum, "compartment": c, "clone": "Other",
                     "count": int(rest.sum()), "freq": float(rest.sum() / counts[c].sum()),
                     "rank": np.nan, "shared": False, "compartments_present": ""})

    sh = freqs[n_present > 1].copy()
    sh.columns = [f"freq_{c}" for c in sh.columns]
    sh.insert(0, "n_compartments", n_present[n_present > 1])
    sh.insert(0, "compartments_present", where[n_present > 1])
    sh.insert(0, "clone", sh.index)
    sh.insert(0, "stratum", stratum)
    sh.insert(0, "patient", patient)
    return summary, flow, sh.reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--concat", required=True, help="concatenated_cdr3_sorted.tsv")
    ap.add_argument("--samplesheet", required=True, help="samplesheet CSV")
    ap.add_argument("--compartment-col", default="origin")
    ap.add_argument("--subject-col", default="patient")
    ap.add_argument("--time-col", default="timepoint",
                    help="cohort time column; per-time strata are skipped if absent")
    ap.add_argument("--fallback-time-col", default="timepoint")
    ap.add_argument("--match-key", default="junction_aa,v_call",
                    help="comma-separated clone key: junction_aa[,v_call]")
    ap.add_argument("--top-n", type=int, default=15, help="clones per compartment in the flow diagram")
    args = ap.parse_args()

    empty = pd.DataFrame()
    meta = pd.read_csv(args.samplesheet, dtype=str).fillna("")
    meta.columns = meta.columns.str.strip()
    write_sample_table(meta, args)

    comp_col = args.compartment_col if args.compartment_col in meta.columns else None
    if comp_col is None:
        write_outputs(empty, empty, empty,
                      f"skipped: samplesheet has no '{args.compartment_col}' column "
                      f"(set --compartment_col to the column naming each sample's compartment)")
        return

    subj_col = next((c for c in (args.subject_col,) + PATIENT_CANDIDATES if c in meta.columns), None)
    if subj_col is None:
        write_outputs(empty, empty, empty, "skipped: samplesheet has no patient/subject column")
        return

    time_col = next((c for c in (args.time_col, args.fallback_time_col)
                     if c and c in meta.columns and (meta[c] != "").any()), None)

    meta = meta[["sample", subj_col, comp_col] + ([time_col] if time_col else [])].copy()
    meta["compartment"] = meta[comp_col].map(normalise_compartment)
    meta = meta[(meta["compartment"] != "") & (meta["compartment"].str.lower() != "unknown")]

    key = [k.strip() for k in args.match_key.split(",") if k.strip()]
    use_v = "v_call" in key
    usecols = ["sample", "junction_aa", "duplicate_count"] + (["v_call"] if use_v else [])
    df = pd.read_csv(args.concat, sep="\t", usecols=usecols, dtype={"sample": str, "junction_aa": str})
    df = df[df["junction_aa"].notna() & (df["junction_aa"] != "")]
    df["v_gene"] = df["v_call"].map(v_gene) if use_v else ""
    df["clone"] = df["junction_aa"] + ("|" + df["v_gene"] if use_v else "")
    df = df.merge(meta, on="sample", how="inner")
    if df.empty:
        write_outputs(empty, empty, empty, "skipped: no clonotype rows matched samplesheet samples with a compartment")
        return

    clone_info = df.drop_duplicates("clone").set_index("clone")[["junction_aa", "v_gene"]]

    summaries, flows, shareds = [], [], []
    for patient, pdf in df.groupby(subj_col):
        strata = [("all", pdf)]
        if time_col:
            strata += [(str(t), tdf) for t, tdf in sorted(pdf.groupby(time_col), key=lambda kv: time_sort_key(kv[0]))
                       if str(t) != ""]
        for stratum, sdf in strata:
            if sdf["compartment"].nunique() < 2:
                continue
            counts = sdf.pivot_table(index="clone", columns="compartment", values="duplicate_count",
                                     aggfunc="sum", fill_value=0)
            n_samples = sdf.groupby("compartment")["sample"].nunique().to_dict()
            s, f, sh = analyse_stratum(counts, str(patient), stratum, args.top_n, n_samples)
            summaries += s
            flows += f
            shareds.append(sh)

    if not summaries:
        write_outputs(empty, empty, empty,
                      f"skipped: no patient was sampled in two or more compartments "
                      f"(compartments found: {', '.join(sorted(meta['compartment'].unique())) or 'none'})")
        return

    flow = pd.DataFrame(flows)
    flow = flow.join(clone_info, on="clone")
    flow.loc[flow["clone"] == "Other", ["junction_aa", "v_gene"]] = ""
    shared = pd.concat(shareds, ignore_index=True).join(clone_info, on="clone")
    write_outputs(pd.DataFrame(summaries), flow, shared,
                  f"ok: {flow['patient'].nunique()} patient(s); compartments "
                  f"{', '.join(sorted(meta['compartment'].unique(), key=compartment_sort_key))}; "
                  f"time column: {time_col or 'none'}")


if __name__ == "__main__":
    sys.exit(main())
