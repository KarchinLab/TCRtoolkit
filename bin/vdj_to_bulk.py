#!/usr/bin/env python3
"""
Bridge 3: VDJ_TO_BULK
Converts VDJ_QC's contigs_after_qc.tsv (per-contig, one row per chain per cell)
into per-sample AIRR-format TSVs for TCRtoolkit — used when there is no GEX object
and TCELL_INTEGRATION is skipped.

Input columns expected (from VDJ_QC notebook):
    sample, barcode, chain, cdr3, v_gene, j_gene
    optionally: patient_id, condition, timepoint, cdr3_nt

Only TRB (beta) chain rows are used.

Usage:
    vdj_to_bulk.py <contigs_after_qc.tsv> [<sample_col>]

Outputs:
    bulk_samples/<sample>_bulk.tsv   — one AIRR TSV per sample
    synthetic_samplesheet.csv        — TCRtoolkit-compatible samplesheet
"""

import sys
import os
import pandas as pd


META_COLS = ['patient_id', 'condition', 'timepoint']


def main():
    if len(sys.argv) < 2:
        print("Usage: vdj_to_bulk.py <contigs_after_qc.tsv> [sample_col]")
        sys.exit(1)

    contigs_file = sys.argv[1]
    sample_col   = sys.argv[2] if len(sys.argv) > 2 else 'sample'

    df = pd.read_csv(contigs_file, sep='\t', low_memory=False)

    required = {'chain', 'cdr3', sample_col}
    missing  = required - set(df.columns)
    if missing:
        raise ValueError(f"contigs_after_qc.tsv is missing columns: {missing}. "
                         f"Available: {list(df.columns)}")

    # keep beta chain only, drop rows with no CDR3
    df = df[df['chain'] == 'TRB'].copy()
    df = df[df['cdr3'].notna() & (df['cdr3'].astype(str) != 'NA')]

    if df.empty:
        raise ValueError("No TRB rows found in contigs_after_qc.tsv after filtering.")

    df['junction_aa'] = df['cdr3'].str.upper()
    df['v_call']      = df['v_gene'] if 'v_gene' in df.columns else None
    df['d_call']      = df['d_gene'] if 'd_gene' in df.columns else None
    df['j_call']      = df['j_gene'] if 'j_gene' in df.columns else None
    df['sequence']    = df['cdr3_nt'] if 'cdr3_nt' in df.columns else df['junction_aa']

    os.makedirs('bulk_samples', exist_ok=True)
    samplesheet_rows = []

    for sample, grp in df.groupby(sample_col):
        agg = (
            grp.groupby(['junction_aa', 'v_call', 'd_call', 'j_call'])
               .agg(
                   duplicate_count=('barcode', 'count'),
                   sequence=('sequence', 'first')
               )
               .reset_index()
        )
        total = agg['duplicate_count'].sum()
        agg['duplicate_frequency_percent'] = (agg['duplicate_count'] / total * 100).round(6)
        agg['sequence_id']       = agg['junction_aa']
        agg['junction']          = agg['sequence']
        agg['junction_aa_length'] = agg['junction_aa'].str.len()
        agg['productive']        = 'true'

        out_cols = [
            'sequence_id', 'junction_aa', 'junction', 'junction_aa_length',
            'v_call', 'd_call', 'j_call',
            'productive', 'duplicate_count', 'duplicate_frequency_percent', 'sequence'
        ]
        out_path = f'bulk_samples/{sample}_bulk.tsv'
        agg[out_cols].to_csv(out_path, sep='\t', index=False)

        row = {'sample': sample, 'file': os.path.abspath(out_path)}
        for col in META_COLS:
            if col in grp.columns:
                val = grp[col].dropna()
                row[col] = val.iloc[0] if len(val) > 0 else ''
        samplesheet_rows.append(row)

    ss = pd.DataFrame(samplesheet_rows)
    ss.to_csv('synthetic_samplesheet.csv', index=False)

    print(f"[Bridge 3] Created {len(samplesheet_rows)} bulk sample files in bulk_samples/")
    print(f"[Bridge 3] Synthetic samplesheet written to synthetic_samplesheet.csv")


if __name__ == '__main__':
    main()
