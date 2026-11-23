#!/usr/bin/env python3
"""
Bridge: SC_REPORT_SAMPLESHEET

The bulk report notebooks (template_qc / discovery_brief / details_*) read a bulk-style
samplesheet: sample, alias, patient, timepoint, origin, file. Single-cell runs have no such
file - their metadata lives in the user's SC sample sheet (keyed by sample) and the
pseudobulk units can be {sample}__{phenotype} rather than plain samples. This builds one
row per pseudobulk unit that passed QC, pulling metadata from the SC sample sheet.

Inputs : units TSV (unit, sample, patient_id, phenotype, file) from the pseudobulk meta,
         optionally the user's SC sample sheet CSV, and optionally the per-cell export
         (full route), whose per-cell metadata fills gaps when the GEX object's sample
         names do not match the sample sheet's.
Output : report_samplesheet.csv
"""

import argparse
import pandas as pd

# Report column -> SC sample-sheet candidates, first match wins.
SOURCES = {
    "alias":     ("alias",),
    "patient":   ("patient", "patient_id", "subject", "subject_id", "donor"),
    "timepoint": ("timepoint", "visit", "day"),
    "origin":    ("origin", "tissue", "compartment", "source"),
    "condition": ("condition", "group", "status"),
}


def main():
    ap = argparse.ArgumentParser(description="SC_REPORT_SAMPLESHEET bridge")
    ap.add_argument("units", help="TSV: unit, sample, patient_id, phenotype, file")
    ap.add_argument("--sample-sheet", default=None, help="user SC sample sheet (CSV)")
    ap.add_argument("--cells", default=None, help="per-cell export TSV (optional)")
    ap.add_argument("--timepoint-col", default="timepoint")
    ap.add_argument("--out", default="report_samplesheet.csv")
    args = ap.parse_args()

    units = pd.read_csv(args.units, sep="\t", dtype=str).fillna("")
    out = pd.DataFrame({"sample": units["unit"], "file": units["file"]})
    base = units["sample"].where(units["sample"] != "", units["unit"])
    out["phenotype"] = units["phenotype"]

    ss = None
    if args.sample_sheet:
        ss = pd.read_csv(args.sample_sheet, dtype=str)
        key = "sample" if "sample" in ss.columns else ss.columns[0]
        ss = ss.drop_duplicates(key).set_index(key)

    sources = dict(SOURCES)
    if args.timepoint_col and args.timepoint_col != "timepoint":
        sources["timepoint"] = (args.timepoint_col,) + SOURCES["timepoint"]

    for col, candidates in sources.items():
        src = next((c for c in candidates if ss is not None and c in ss.columns), None)
        out[col] = base.map(ss[src]).fillna("") if src else ""

    # Gap-fill from the per-cell export: most common value per sample.
    if args.cells:
        cells = pd.read_csv(args.cells, sep="\t", dtype=str, low_memory=False)
        if "sample" in cells.columns:
            for col, candidates in sources.items():
                src = next((c for c in candidates if c in cells.columns), None)
                if not src or col == "alias":
                    continue
                per_sample = (cells[["sample", src]].dropna()
                              .loc[lambda d: d[src].str.strip() != ""]
                              .groupby("sample")[src].agg(lambda v: v.mode().iloc[0]))
                fill = base.map(per_sample).fillna("")
                out[col] = out[col].where(out[col] != "", fill)

    # Fallbacks the notebooks need non-empty: patient groups per-subject tabs, alias labels
    # every plot, origin splits tissue/blood panels. Phenotype units keep their suffix in
    # the alias so the per-phenotype views stay distinguishable.
    out["patient"] = out["patient"].where(out["patient"] != "", units["patient_id"])
    out["patient"] = out["patient"].where(out["patient"] != "", base)
    alias = out["alias"].where(out["alias"] != "", base)
    out["alias"] = alias.where(out["phenotype"] == "", alias + "__" + out["phenotype"])
    out["origin"] = out["origin"].where(out["origin"] != "", "unknown")

    # Extra user columns ride along (they become groupable metadata in the notebooks).
    if ss is not None:
        skip = {"path", "file"} | {c for cands in sources.values() for c in cands}
        for c in ss.columns:
            if c not in skip and c not in out.columns:
                out[c] = base.map(ss[c]).fillna("")

    lead = ["sample", "alias", "patient", args.timepoint_col or "timepoint", "origin", "file"]
    if args.timepoint_col and args.timepoint_col != "timepoint":
        out = out.rename(columns={"timepoint": args.timepoint_col})
    out = out[lead + [c for c in out.columns if c not in lead]]
    out.sort_values("sample").to_csv(args.out, index=False)
    print(f"[SC_REPORT_SAMPLESHEET] {len(out)} unit(s) -> {args.out}")


if __name__ == "__main__":
    main()
