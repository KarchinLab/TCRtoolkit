#!/usr/bin/env python3
"""
Bridge: SC_TO_CDR3
Converts TCELL_INTEGRATION per-cell export_cells.tsv into the canonical TCRtoolkit
AIRR clonotype schema (per pseudobulk unit), bypassing ANNOTATE_PROCESS.

Pooling policy (integration design decision):
  * Default — pool by SAMPLE. Each sample is one unit; clonotypes are grouped by
    (junction_aa, v_call, j_call). Patient is carried as metadata so the shared
    PATIENT step pools per-patient for GIANA/GLIPH2/tcrdist clustering, and
    phenotype (annot) enters downstream via CLUSTER_TO_SC / CONGA / repertoire.
  * Optional (--by-phenotype) — secondary per-cell-type view. Units become
    (sample x annot), id = "{sample}__{annot}". Under-powered units are dropped
    later by PSEUDOBULK_QC (min clones/cells). Off by default.

Schema conformance: output columns + ORDER match ANNOTATE_PROCESS so the shared
bulk engine runs unmodified (positional sort/dedup: $1=junction_aa, $2=v_call):

    junction_aa, v_call, d_call, j_call, duplicate_count,
    junction_aa_length, duplicate_frequency_percent,
    sequence, sequence_id, junction, sample

Outputs:
    per_sample/{unit}_cdr3.tsv  — one file per pseudobulk unit
    concat_cdr3.tsv             — all units concatenated (input for ANNOTATE_SORT_CDR3)
    unit_map.csv                — unit -> sample/patient/phenotype/file (for meta building)
"""

import sys
import os
import argparse
import pandas as pd

OUTPUT_COLS = [
    "junction_aa", "v_call", "d_call", "j_call",
    "duplicate_count", "junction_aa_length", "duplicate_frequency_percent",
    "sequence", "sequence_id", "junction", "sample",
]


def main():
    ap = argparse.ArgumentParser(description="SC_TO_CDR3 bridge")
    ap.add_argument("export_cells", help="tcr_export_cells_with_embedding.tsv")
    ap.add_argument("--by-phenotype", action="store_true",
                    help="stratify pseudobulk units by cell-type (annot); off by default")
    ap.add_argument("--pheno-col", default="annot", help="phenotype column name")
    args = ap.parse_args()

    df = pd.read_csv(args.export_cells, sep='\t', low_memory=False)

    # Resolve CDR3b amino-acid sequence → junction_aa
    if 'cdr3b' in df.columns:
        df['junction_aa'] = df['cdr3b']
    elif 'junction_aa' in df.columns:
        df['junction_aa'] = df['junction_aa']
    else:
        raise ValueError("Export file must contain 'cdr3b' or 'junction_aa'.")

    # Resolve V / J gene calls, ensuring an IMGT allele suffix (e.g. TRBV10-3 -> TRBV10-3*01).
    # tcrdist3's reference db is keyed by allele and cannot map CDR1/CDR2 without one.
    def add_allele(g):
        if g is None or pd.isna(g) or str(g) in ('', 'None'):
            return g
        s = str(g)
        return s if '*' in s else f"{s}*01"

    df['v_call'] = (df['trbv'].apply(add_allele) if 'trbv' in df.columns else None)
    df['j_call'] = (df['trbj'].apply(add_allele) if 'trbj' in df.columns else None)

    if 'sample' not in df.columns:
        raise ValueError("Export file must contain a 'sample' column.")

    # Patient metadata (carried so PATIENT can pool per-patient); default to sample if absent.
    df['patient'] = df['patient'] if 'patient' in df.columns else df['sample']

    # Phenotype stratification (optional)
    by_pheno = args.by_phenotype and args.pheno_col in df.columns
    if args.by_phenotype and not by_pheno:
        print(f"[SC_TO_CDR3] --by-phenotype requested but column '{args.pheno_col}' "
              f"not found; falling back to sample-level pooling.")
    df['phenotype'] = df[args.pheno_col].astype(str) if by_pheno else ''
    df['unit'] = (df['sample'].astype(str) + '__' + df['phenotype']) if by_pheno else df['sample'].astype(str)

    # Keep only cells with a valid beta CDR3
    df = df[df['junction_aa'].notna() & (df['junction_aa'] != '') & (df['junction_aa'] != 'None')]
    print(f"[SC_TO_CDR3] {len(df)} cells with valid junction_aa across "
          f"{df['unit'].nunique()} unit(s) (by_phenotype={by_pheno}).")

    os.makedirs('per_sample', exist_ok=True)
    all_frames = []
    unit_rows = []

    for unit, grp in df.groupby('unit'):
        agg = (
            grp.groupby(['junction_aa', 'v_call', 'j_call'], dropna=False)
               .size()
               .reset_index(name='duplicate_count')
        )
        agg['d_call']      = ''
        total = agg['duplicate_count'].sum()
        agg['duplicate_frequency_percent'] = (agg['duplicate_count'] / total * 100).round(6) if total else 0.0
        agg['junction_aa_length'] = agg['junction_aa'].astype(str).str.len()
        agg['sequence']    = ''
        agg['sequence_id'] = agg['junction_aa']
        agg['junction']    = ''
        agg['sample']      = unit
        agg = agg[OUTPUT_COLS]

        out_path = f'per_sample/{unit}_cdr3.tsv'
        agg.to_csv(out_path, sep='\t', index=False)
        all_frames.append(agg)

        # First non-null patient / phenotype for this unit
        patient = grp['patient'].dropna().astype(str)
        patient = patient.iloc[0] if len(patient) else ''
        phenotype = grp['phenotype'].iloc[0] if len(grp) else ''
        unit_rows.append({'sample': unit, 'patient': patient,
                          'phenotype': phenotype, 'file': os.path.abspath(out_path)})
        print(f"[SC_TO_CDR3]   {unit}: {len(agg)} clonotypes (patient={patient}).")

    if not all_frames:
        raise ValueError("[SC_TO_CDR3] No clonotypes found after filtering.")

    concat = pd.concat(all_frames, ignore_index=True)
    concat.to_csv('concat_cdr3.tsv', sep='\t', index=False)
    pd.DataFrame(unit_rows).to_csv('unit_map.csv', index=False)
    print(f"[SC_TO_CDR3] concat_cdr3.tsv ({len(concat)} rows), unit_map.csv ({len(unit_rows)} units) written.")


if __name__ == '__main__':
    main()
