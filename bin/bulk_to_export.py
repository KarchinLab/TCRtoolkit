#!/usr/bin/env python3
"""
Bridge: BULK_TO_EXPORT

Builds a per-cell export table (the format REPERTOIRE / MASTER_SUMMARY consume) from
pseudobulk clonotype data, so those reports can run in the VDJ-only route where there is
no GEX Seurat object.

The repertoire report derives per-clone cell counts by COUNTING ROWS grouped by
(sample, clone_id) — i.e. it treats each row as one cell. So each clonotype with
duplicate_count = N is expanded into N rows (one synthetic "cell" each). Without a GEX
object there is no cell-type annotation, so `annot` is a constant ("Unannotated").

Input : concatenated canonical clonotype table (junction_aa, v_call, j_call,
        duplicate_count, ..., sample) — e.g. ANNOTATE_FROM_CONCAT's concat_cdr3_sorted.
Optional: samplesheet CSV mapping sample -> patient_id (else patient_id = sample).

Output: export_cells.tsv with columns REPERTOIRE/MASTER_SUMMARY resolve:
        cell_id, sample, patient_id, clone_id, clone_size, annot, has_tcr, paired_tcr
"""

import sys
import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description="BULK_TO_EXPORT bridge")
    ap.add_argument("concat_cdr3", help="concatenated canonical clonotype TSV (has 'sample')")
    ap.add_argument("--samplesheet", default=None, help="CSV with sample,patient_id (optional)")
    ap.add_argument("--out", default="export_cells.tsv")
    args = ap.parse_args()

    df = pd.read_csv(args.concat_cdr3, sep="\t", low_memory=False)

    for col in ("junction_aa", "duplicate_count", "sample"):
        if col not in df.columns:
            sys.exit(f"[BULK_TO_EXPORT] input missing required column '{col}'. Have: {list(df.columns)}")

    # Clone identity = CDR3b + V gene (falls back to junction_aa alone if v_call absent)
    if "v_call" in df.columns:
        df["clone_id"] = df["junction_aa"].astype(str) + "_" + df["v_call"].astype(str)
    else:
        df["clone_id"] = df["junction_aa"].astype(str)

    df["clone_size"] = pd.to_numeric(df["duplicate_count"], errors="coerce").fillna(0).astype(int)
    df = df[df["clone_size"] > 0].copy()
    if df.empty:
        sys.exit("[BULK_TO_EXPORT] no clonotypes with positive counts.")

    # patient_id: from samplesheet if given, else = sample
    if args.samplesheet:
        ss = pd.read_csv(args.samplesheet)
        pcol = next((c for c in ("patient_id", "patient") if c in ss.columns), None)
        if pcol and "sample" in ss.columns:
            pmap = dict(zip(ss["sample"].astype(str), ss[pcol].astype(str)))
            df["patient_id"] = df["sample"].astype(str).map(pmap).fillna(df["sample"].astype(str))
        else:
            df["patient_id"] = df["sample"].astype(str)
    else:
        df["patient_id"] = df["sample"].astype(str)

    # Expand each clonotype into `clone_size` per-cell rows (one row = one cell)
    expanded = df.loc[df.index.repeat(df["clone_size"])].reset_index(drop=True)
    expanded["annot"]      = "Unannotated"   # no GEX → no cell-type label
    expanded["has_tcr"]    = "TRUE"
    expanded["paired_tcr"] = "FALSE"
    expanded["cell_id"]    = ["cell_%d" % i for i in range(len(expanded))]

    out_cols = ["cell_id", "sample", "patient_id", "clone_id", "clone_size",
                "annot", "has_tcr", "paired_tcr"]
    expanded[out_cols].to_csv(args.out, sep="\t", index=False)
    print(f"[BULK_TO_EXPORT] {df['clone_id'].nunique()} clonotypes across "
          f"{df['sample'].nunique()} sample(s) -> {len(expanded)} cell rows in {args.out}")


if __name__ == "__main__":
    main()
