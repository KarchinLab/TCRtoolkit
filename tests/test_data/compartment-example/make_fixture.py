#!/usr/bin/env python3
"""
Builds the compartment-example test fixture from 10x Cell Ranger AIRR rearrangement
tables (GBM CAN-3110 trial, Ling et al. 2025): two patients sampled in tumor and blood at
matched protocol days.

Each 10x table (one row per contig) is collapsed to a clonotype-level AIRR table: productive
TRB contigs, one row per (junction_aa, v_call, j_call), duplicate_count = number of cells.
Clones are down-sampled (the largest kept, plus a random set of the rest) so the fixture
stays small while keeping clones shared between compartments.

Usage: make_fixture.py <AIRR_files dir> <out dir>
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SAMPLES = [
    # sample, file stem, patient, origin, timepoint (binned cohort day), protocol_day (actual day)
    ("P002_d0_Tumor",  "GBM1-DFCI2_S1-TCR-LIB",           "P002", "Tumor", "0",  "d0"),
    ("P002_d30_Tumor", "GBM1-DFCI2_S3-TCR-LIB",           "P002", "Tumor", "30", "d30"),
    ("P002_d0_Blood",  "BTC-GBM1-DFCI2-S1-71-S1-TCR-LIB", "P002", "PBMC", "0",  "d0"),
    ("P002_d30_Blood", "BTC-GBM1-DFCI2-S3-71-RK-TCR-LIB", "P002", "PBMC", "30", "d30"),
    ("P004_d0_Tumor",  "GBM1-DFCI4-S1-TCR-LIB",           "P004", "Tumor", "0",  "d0"),
    ("P004_d33_Tumor", "GBM1-DFCI4-S3-L1-TCR-LIB",        "P004", "Tumor", "30", "d33"),
    ("P004_d0_Blood",  "BTC-GBM1-DFCI4-S1-73-S1-TCR-LIB", "P004", "PBMC", "0",  "d0"),
    ("P004_d33_Blood", "BTC-GBM1-DFCI4-S3-73-S3-TCR-LIB", "P004", "PBMC", "30", "d33"),
]
TOP, RANDOM = 150, 250
COLS = ["sequence_id", "sequence", "junction", "junction_aa", "v_call", "d_call", "j_call",
        "duplicate_count", "junction_aa_length", "duplicate_frequency_percent", "productive"]


def collapse(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str, low_memory=False)
    df = df[(df["productive"].str.upper().isin(["T", "TRUE"])) & df["v_call"].str.startswith("TRB", na=False)
            & df["junction_aa"].notna()]
    # First call of an ambiguous list; Cell Ranger omits the allele, which AIRR consumers
    # (e.g. tcrdist3's gene lookup) expect, so add *01.
    for c in ("v_call", "d_call", "j_call"):
        g = df[c].fillna("").str.split(",").str[0]
        df[c] = g.where((g == "") | g.str.contains(r"\*", regex=True), g + "*01")
    g = (df.groupby(["junction_aa", "v_call", "j_call"])
           .agg(duplicate_count=("cell_id", "nunique"), d_call=("d_call", "first"),
                junction=("junction", "first"), sequence=("sequence", "first"))
           .reset_index())
    return g


def main(src: str, out: str):
    rng = np.random.default_rng(7)
    out = Path(out)
    tables = {s: collapse(Path(src) / f"{stem}_airr_rearrangement.tsv") for s, stem, *_ in SAMPLES}
    # Keep every clone present in both compartments of a patient so the flows are not empty.
    for patient in {p for _, _, p, *_ in SAMPLES}:
        names = [s for s, _, p, *_ in SAMPLES if p == patient]
        by_origin = {}
        for s, _, p, o, *_ in SAMPLES:
            if p == patient:
                by_origin.setdefault(o, set()).update(tables[s]["junction_aa"])
        shared = set.intersection(*by_origin.values())
        for s in names:
            t = tables[s].sort_values("duplicate_count", ascending=False)
            keep = t.head(TOP)
            rest = t.iloc[TOP:]
            keep = pd.concat([keep, rest[rest["junction_aa"].isin(shared)],
                              rest.sample(min(RANDOM, len(rest)), random_state=7)]).drop_duplicates()
            tables[s] = keep
    rows = []
    for s, stem, patient, origin, tp, day in SAMPLES:
        t = tables[s].copy()
        t["sequence_id"] = [f"{s}_{i}" for i in range(len(t))]
        t["junction_aa_length"] = t["junction_aa"].str.len()
        t["duplicate_frequency_percent"] = 100 * t["duplicate_count"] / t["duplicate_count"].sum()
        t["productive"] = True
        t[COLS].to_csv(out / f"{s}.tsv", sep="\t", index=False)
        rows.append({"sample": s, "alias": s, "patient": patient, "timepoint": tp,
                     "protocol_day": day, "origin": origin,
                     # Relative to the repository root, like params.yml (run from there).
                     "file": f"tests/test_data/compartment-example/{s}.tsv"})
        print(s, len(t), "clones")
    pd.DataFrame(rows).to_csv(out / "samplesheet.csv", index=False)


if __name__ == "__main__":
    main(*sys.argv[1:3])
