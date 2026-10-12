#!/usr/bin/env python3
"""
Bridge: BULK_TO_EXPORT

Builds a per-cell export table (the format CLUSTER_TO_SC / CONSENSUS / REPERTOIRE /
MASTER_SUMMARY consume) from clonotype-level data, so those steps can run without a GEX
Seurat object - on the single-cell VDJ-only route and on bulk TCR data.

Two granularities (--level):
  cell       Each clonotype with duplicate_count = N is expanded into N rows (one synthetic
             "cell" each). The repertoire report derives per-clone cell counts by COUNTING
             ROWS grouped by (sample, clone_id), so this is what the single-cell VDJ-only
             route uses: there, duplicate_count really is a cell count.
  clonotype  One row per (sample, clonotype), clone_size = summed duplicate_count. Used for
             bulk data, where duplicate_count is a read/template count that can run into the
             millions per sample, so per-count expansion is not feasible. Downstream "cell"
             counts then mean clonotypes.

Without a GEX object there is no cell-type annotation, so `annot` is a constant
("Unannotated").

Input : concatenated canonical clonotype table (junction_aa, v_call, j_call,
        duplicate_count, ..., sample) - e.g. ANNOTATE's concat_cdr3_sorted.
Optional: samplesheet CSV mapping sample -> patient / timepoint / origin / condition.

Output: export_cells.tsv with columns
        cell_id, sample, patient_id, timepoint, condition, origin, clone_id, clone_size,
        junction_aa, v_call, j_call, annot, has_tcr, paired_tcr
"""

import sys
import argparse
import pandas as pd

# Samplesheet column -> export column, first match wins.
META_SOURCES = {
    "patient_id": ("patient_id", "patient", "subject_id", "subject"),
    "timepoint":  ("timepoint", "visit", "day"),
    "condition":  ("condition", "group", "status"),
    "origin":     ("origin", "tissue", "compartment", "source"),
}


def main():
    ap = argparse.ArgumentParser(description="BULK_TO_EXPORT bridge")
    ap.add_argument("concat_cdr3", help="concatenated canonical clonotype TSV (has 'sample')")
    ap.add_argument("--samplesheet", default=None, help="CSV with sample + metadata columns (optional)")
    ap.add_argument("--level", choices=("cell", "clonotype"), default="cell")
    ap.add_argument("--out", default="export_cells.tsv")
    args = ap.parse_args()

    df = pd.read_csv(args.concat_cdr3, sep="\t", low_memory=False)

    for col in ("junction_aa", "duplicate_count", "sample"):
        if col not in df.columns:
            sys.exit(f"[BULK_TO_EXPORT] input missing required column '{col}'. Have: {list(df.columns)}")

    df["sample"] = df["sample"].astype(str)
    for col in ("v_call", "j_call"):
        df[col] = df[col].fillna("").astype(str) if col in df.columns else ""

    # Clone identity = CDR3b + V gene (falls back to junction_aa alone if v_call absent)
    df["clone_id"] = df["junction_aa"].astype(str)
    has_v = df["v_call"].ne("")
    df.loc[has_v, "clone_id"] = df.loc[has_v, "junction_aa"].astype(str) + "_" + df.loc[has_v, "v_call"]

    df["clone_size"] = pd.to_numeric(df["duplicate_count"], errors="coerce").fillna(0).astype(int)
    df = df[df["clone_size"] > 0].copy()
    if df.empty:
        sys.exit("[BULK_TO_EXPORT] no clonotypes with positive counts.")

    # A clonotype can appear on several rows of one sample (different nucleotide sequences
    # encoding the same CDR3aa + V); collapse them so each (sample, clone_id) is one unit.
    df = (df.groupby(["sample", "clone_id"], sort=False, as_index=False)
            .agg(clone_size=("clone_size", "sum"), junction_aa=("junction_aa", "first"),
                 v_call=("v_call", "first"), j_call=("j_call", "first")))

    # Per-sample metadata from the samplesheet; patient_id falls back to sample.
    meta = pd.DataFrame({"sample": df["sample"].unique()})
    if args.samplesheet:
        ss = pd.read_csv(args.samplesheet, dtype=str)
        if "sample" in ss.columns:
            ss = ss.drop_duplicates("sample")
            for out_col, candidates in META_SOURCES.items():
                src = next((c for c in candidates if c in ss.columns), None)
                if src:
                    meta[out_col] = meta["sample"].map(dict(zip(ss["sample"], ss[src])))
    for out_col in META_SOURCES:
        if out_col not in meta.columns:
            meta[out_col] = pd.NA
    meta["patient_id"] = meta["patient_id"].fillna(meta["sample"])
    df = df.merge(meta, on="sample", how="left")

    if args.level == "cell":
        # Expand each clonotype into `clone_size` per-cell rows (one row = one cell)
        out = df.loc[df.index.repeat(df["clone_size"])].reset_index(drop=True)
        out["cell_id"] = ["cell_%d" % i for i in range(len(out))]
    else:
        out = df.reset_index(drop=True)
        out["cell_id"] = ["clonotype_%d" % i for i in range(len(out))]

    out["annot"]      = "Unannotated"   # no GEX -> no cell-type label
    out["has_tcr"]    = "TRUE"
    out["paired_tcr"] = "FALSE"

    out_cols = ["cell_id", "sample", "patient_id", "timepoint", "condition", "origin",
                "clone_id", "clone_size", "junction_aa", "v_call", "j_call",
                "annot", "has_tcr", "paired_tcr"]
    out[out_cols].to_csv(args.out, sep="\t", index=False)
    unit = "cell" if args.level == "cell" else "clonotype"
    print(f"[BULK_TO_EXPORT] {df['clone_id'].nunique()} clonotypes across "
          f"{df['sample'].nunique()} sample(s) -> {len(out)} {unit} rows in {args.out}")


if __name__ == "__main__":
    main()
