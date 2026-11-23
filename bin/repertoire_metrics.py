#!/usr/bin/env python3
"""
REPERTOIRE_METRICS

The one place repertoire-level metrics are computed, for both modalities. The report
notebooks and the Master Summary read these tables instead of each recomputing them.

Clone unit:
    bulk         a TRB clonotype row of the concatenated CDR3 table (junction_aa + V/D/J)
                 for diversity; its junction_aa for overlap and sharing. These are the
                 definitions the bulk engine has always used.
    single-cell  a paired clone (clone_id 'A:<CDR3a>|B:<CDR3b>', first chain of each
                 type per cell), counted in cells. Alpha- or beta-only cells keep their
                 single-chain id.

Inputs:
    --clones        bulk: concatenated_cdr3_sorted.tsv (sample, junction_aa, v_call,
                    d_call, j_call, duplicate_count)
                    single-cell: paired_clones.tsv (sample, clone_id, n_cells)
    --sample-stats  SAMPLE_CALC's per-sample stats (productive, CDR3 length, convergence);
                    the diversity columns are added here. Its samples define the cohort.
    --samplesheet   per-sample metadata (patient, timepoint, condition, origin)

Outputs (--outdir):
    sample_stats.csv        sample stats + diversity (num_clones, num_TCRs, clonality,
                            simpson_index[_corrected], shannon_entropy, inverse_simpson,
                            gini_coefficient, hill_q0..q2) + clone_unit
    diversity_by_group.tsv  the same metrics on clones pooled per patient / condition /
                            timepoint / origin
    clone_size_bins.tsv     per sample x frequency bin: clones and TCRs (cells or templates)
    overlap_samples.tsv     every sample pair: shared clones, Jaccard, Sorensen, Morisita-Horn
    overlap_patients.tsv    every patient pair (clones pooled per patient): shared, Jaccard
    clone_sharing.tsv       per clone: samples / patients / conditions it occurs in
    repertoire_summary.tsv  metric / value rollup
"""

import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import sparse

GROUP_LEVELS = ["patient", "condition", "timepoint", "origin"]
# Sample-sheet column aliases, first match wins.
META_ALIASES = {
    "patient":   ["patient", "patient_id", "subject_id"],
    "condition": ["condition"],
    "timepoint": ["timepoint"],
    "origin":    ["origin", "tissue"],
}
DIVERSITY_COLS = ["num_clones", "num_TCRs", "simpson_index", "simpson_index_corrected",
                  "clonality", "shannon_entropy", "inverse_simpson", "gini_coefficient",
                  "hill_q0", "hill_q1", "hill_q2"]


def diversity(counts):
    """Diversity of one repertoire from its clone counts (definitions as in sample_calc)."""
    c = np.asarray(counts, dtype=float)
    c = c[c > 0]
    n, total = len(c), c.sum()
    if n == 0:
        return {k: np.nan for k in DIVERSITY_COLS}
    p = c / total
    h_nat = -np.sum(p * np.log(p))
    h_2 = h_nat / np.log(2)
    sum_p2 = np.sum(p ** 2)
    sorted_p = np.sort(p)
    idx = np.arange(1, n + 1)
    return {
        "num_clones": int(n),
        "num_TCRs": int(total) if float(total).is_integer() else total,
        "simpson_index": sum_p2,
        "simpson_index_corrected": (np.sum(c * (c - 1)) / (total * (total - 1))
                                    if total > 1 else np.nan),
        "clonality": 1 - h_2 / np.log2(n) if n > 1 else np.nan,
        "shannon_entropy": h_nat,
        "inverse_simpson": 1.0 / sum_p2,
        "gini_coefficient": ((2 * idx - n - 1) * sorted_p).sum() / (n * sorted_p.sum()),
        "hill_q0": int(n),
        "hill_q1": np.exp(h_nat),
        "hill_q2": 1.0 / sum_p2,
    }


def load_clones(path):
    """Return (long table sample/div_key/overlap_key/count, clone_unit)."""
    df = pd.read_csv(path, sep="\t", low_memory=False)
    if "clone_id" in df.columns:
        df = df.dropna(subset=["sample", "clone_id"])
        df = df[df["clone_id"].astype(str).str.strip() != ""]
        count = df["n_cells"] if "n_cells" in df.columns else pd.Series(1, index=df.index)
        out = pd.DataFrame({"sample": df["sample"].astype(str),
                            "div_key": df["clone_id"].astype(str),
                            "count": pd.to_numeric(count, errors="coerce").fillna(0)})
        out["overlap_key"] = out["div_key"]
        unit = "paired_clone"
    elif "junction_aa" in df.columns:
        df = df.dropna(subset=["sample", "junction_aa"])
        parts = [df[c].fillna("").astype(str) for c in ("junction_aa", "v_call", "d_call", "j_call")
                 if c in df.columns]
        out = pd.DataFrame({"sample": df["sample"].astype(str),
                            "div_key": parts[0].str.cat(parts[1:], sep="|") if len(parts) > 1 else parts[0],
                            "overlap_key": df["junction_aa"].astype(str),
                            "count": pd.to_numeric(df["duplicate_count"], errors="coerce").fillna(0)})
        unit = "trb_clonotype"
    else:
        sys.exit("[repertoire_metrics] --clones needs a clone_id or junction_aa column")
    # Rows are kept as they are: a bulk row is one clonotype even when another row of the
    # sample shares its amino-acid CDR3 and genes (different nucleotides), as sample_calc
    # has always counted them. Pooled and pairwise metrics aggregate by key themselves.
    return out[out["count"] > 0].reset_index(drop=True), unit


def load_meta(path, samples):
    meta = pd.DataFrame({"sample": samples})
    if path and os.path.isfile(path) and os.path.getsize(path) > 0:
        ss = pd.read_csv(path, sep=None, engine="python", dtype=str, encoding="utf-8-sig")
        if "sample" in ss.columns:
            keep = {"sample": "sample"}
            for level, aliases in META_ALIASES.items():
                col = next((a for a in aliases if a in ss.columns), None)
                if col:
                    keep[col] = level
            ss = ss[list(keep)].rename(columns=keep).drop_duplicates("sample")
            meta = meta.merge(ss, on="sample", how="left")
    for level in GROUP_LEVELS:
        if level not in meta.columns:
            meta[level] = np.nan
        meta[level] = meta[level].replace("", np.nan)
    return meta


def pairwise(clones, unit_col, key="overlap_key"):
    """Shared / Jaccard / Sorensen / Morisita-Horn for every pair of units (long table)."""
    agg = clones.groupby([unit_col, key], as_index=False)["count"].sum()
    units = sorted(agg[unit_col].unique())
    keys = pd.Index(agg[key].unique())
    ui = {u: i for i, u in enumerate(units)}
    M = sparse.csr_matrix((agg["count"].to_numpy(float),
                           (agg[unit_col].map(ui).to_numpy(), keys.get_indexer(agg[key]))),
                          shape=(len(units), len(keys)))
    B = M.copy()
    B.data[:] = 1
    shared = (B @ B.T).toarray()
    size = np.asarray(B.sum(axis=1)).ravel()
    union = size[:, None] + size[None, :] - shared
    jac = np.divide(shared, union, out=np.zeros_like(shared, dtype=float), where=union > 0)
    sor = np.divide(2 * shared, size[:, None] + size[None, :],
                    out=np.zeros_like(shared, dtype=float), where=(size[:, None] + size[None, :]) > 0)
    tot = np.asarray(M.sum(axis=1)).ravel()
    lam = np.divide(np.asarray(M.multiply(M).sum(axis=1)).ravel(), tot ** 2,
                    out=np.zeros_like(tot), where=tot > 0)
    prod = (M @ M.T).toarray()
    den = (lam[:, None] + lam[None, :]) * np.outer(tot, tot)
    mh = np.divide(2 * prod, den, out=np.zeros_like(prod), where=den > 0)
    a, b = np.meshgrid(np.arange(len(units)), np.arange(len(units)), indexing="ij")
    a, b = a.ravel(), b.ravel()
    return pd.DataFrame({"a": np.asarray(units)[a], "b": np.asarray(units)[b],
                         "n_shared": shared[a, b].astype(int), "jaccard": jac[a, b],
                         "sorensen": sor[a, b], "morisita_horn": mh[a, b]})


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--clones", required=True)
    ap.add_argument("--sample-stats", default=None)
    ap.add_argument("--samplesheet", default=None)
    ap.add_argument("--bin-breaks", default="0;1e-4;0.001;0.01;0.1;1.0")
    ap.add_argument("--bin-labels", default="Rare;Small;Medium;Large;Hyperexpanded")
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    def out(name):
        return os.path.join(args.outdir, name)

    clones, unit = load_clones(args.clones)

    stats = None
    if args.sample_stats and os.path.isfile(args.sample_stats) and os.path.getsize(args.sample_stats) > 0:
        stats = pd.read_csv(args.sample_stats)
        stats["sample"] = stats["sample"].astype(str)
        # Diversity is owned here; drop any stale copy so there is a single definition.
        stats = stats.drop(columns=[c for c in DIVERSITY_COLS + ["clone_unit"] if c in stats.columns])
        clones = clones[clones["sample"].isin(set(stats["sample"]))]
    samples = sorted(clones["sample"].unique())
    if not samples:
        sys.exit("[repertoire_metrics] no clones left for the samples in --sample-stats")
    meta = load_meta(args.samplesheet, samples)
    clones = clones.merge(meta, on="sample", how="left")
    print(f"[repertoire_metrics] {len(samples)} samples, {len(clones)} sample-clone rows, unit={unit}")

    # ── Per-sample diversity → sample_stats.csv ─────────────────────────────
    div = pd.DataFrame([{"sample": s, **diversity(g["count"])}
                        for s, g in clones.groupby("sample")])
    div["clone_unit"] = unit
    if stats is not None:
        rest = [c for c in stats.columns if c != "sample"]
        merged = stats.merge(div, on="sample", how="left")
        merged = merged[["sample"] + DIVERSITY_COLS[:5] + rest + DIVERSITY_COLS[5:] + ["clone_unit"]]
    else:
        merged = div
    merged.sort_values("sample").to_csv(out("sample_stats.csv"), index=False)

    # ── Pooled diversity per group ──────────────────────────────────────────
    rows = []
    for level in GROUP_LEVELS:
        sub = clones.dropna(subset=[level])
        for grp, g in sub.groupby(level):
            pooled = g.groupby("div_key")["count"].sum()
            rows.append({"level": level, "group": grp, "n_samples": g["sample"].nunique(),
                         **diversity(pooled)})
    pd.DataFrame(rows, columns=["level", "group", "n_samples"] + DIVERSITY_COLS) \
        .assign(clone_unit=unit).to_csv(out("diversity_by_group.tsv"), sep="\t", index=False)

    # ── Clone-size bins ─────────────────────────────────────────────────────
    breaks = [float(x) for x in args.bin_breaks.split(";")]
    labels = args.bin_labels.split(";")
    c = clones[["sample", "div_key", "count"]].copy()
    c["freq"] = c["count"] / c.groupby("sample")["count"].transform("sum")
    c["clone_bin"] = pd.cut(c["freq"], bins=breaks, labels=labels, include_lowest=True, right=True)
    bins = (c.groupby(["sample", "clone_bin"], observed=False)
              .agg(n_clones=("div_key", "size"), n_tcrs=("count", "sum")).reset_index())
    bins["frac_clones"] = bins["n_clones"] / bins.groupby("sample")["n_clones"].transform("sum")
    bins["frac_tcrs"] = bins["n_tcrs"] / bins.groupby("sample")["n_tcrs"].transform("sum")
    bins["clone_unit"] = unit
    bins.to_csv(out("clone_size_bins.tsv"), sep="\t", index=False)

    # ── Overlap ─────────────────────────────────────────────────────────────
    ov = pairwise(clones, "sample").rename(columns={"a": "sample_a", "b": "sample_b"})
    pmap = meta.set_index("sample")["patient"].to_dict()
    ov.insert(2, "patient_a", ov["sample_a"].map(pmap))
    ov.insert(3, "patient_b", ov["sample_b"].map(pmap))
    ov["clone_unit"] = unit
    ov.to_csv(out("overlap_samples.tsv"), sep="\t", index=False)

    pc = clones.dropna(subset=["patient"])
    if pc["patient"].nunique() > 0:
        # Patient overlap is on clone presence, as the reports have always shown it.
        po = pairwise(pc.assign(count=1), "patient")[["a", "b", "n_shared", "jaccard"]]
        po.rename(columns={"a": "patient_a", "b": "patient_b"}).assign(clone_unit=unit) \
          .to_csv(out("overlap_patients.tsv"), sep="\t", index=False)

    # ── Sharing ─────────────────────────────────────────────────────────────
    share = (clones.groupby("overlap_key")
                   .agg(n_samples=("sample", "nunique"), n_patients=("patient", "nunique"),
                        n_conditions=("condition", "nunique"), n_timepoints=("timepoint", "nunique"),
                        total_count=("count", "sum"),
                        samples=("sample", lambda s: ",".join(sorted(set(s)))))
                   .reset_index().rename(columns={"overlap_key": "clone"})
                   .sort_values(["n_samples", "n_patients", "total_count"], ascending=False))
    share["clone_unit"] = unit
    share.to_csv(out("clone_sharing.tsv"), sep="\t", index=False)

    # ── Rollup ──────────────────────────────────────────────────────────────
    summary = [
        ("Clone unit", unit),
        ("Samples", len(samples)),
        ("Patients", meta["patient"].nunique()),
        ("Conditions", meta["condition"].nunique()),
        ("Timepoints", meta["timepoint"].nunique()),
        ("Unique clones", clones["overlap_key"].nunique()),
        ("Total TCRs (cells or templates)", int(clones["count"].sum())),
        ("Clones in >1 sample", int((share["n_samples"] > 1).sum())),
        ("Clones in >1 patient", int((share["n_patients"] > 1).sum())),
        ("Median clonality", round(float(div["clonality"].median()), 4)),
    ]
    pd.DataFrame(summary, columns=["metric", "value"]) \
      .to_csv(out("repertoire_summary.tsv"), sep="\t", index=False)
    print("[repertoire_metrics] done")


if __name__ == "__main__":
    main()
