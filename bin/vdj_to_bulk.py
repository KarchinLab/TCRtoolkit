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
    import argparse
    ap = argparse.ArgumentParser(description="VDJ_TO_BULK bridge")
    ap.add_argument("contigs_after_qc")
    ap.add_argument("sample_col", nargs="?", default="sample")
    ap.add_argument("--sample-sheet", default=None,
                    help="user SC sample sheet (CSV); patient_id/condition/timepoint are "
                         "joined from here by sample when the contigs lack them")
    args = ap.parse_args()

    contigs_file = args.contigs_after_qc
    sample_col   = args.sample_col

    # Optional per-sample metadata from the user's sample sheet (keyed by 'sample').
    meta_lookup = {}
    if args.sample_sheet:
        ss_in = pd.read_csv(args.sample_sheet)
        key = 'sample' if 'sample' in ss_in.columns else ss_in.columns[0]
        for _, r in ss_in.iterrows():
            meta_lookup[str(r[key])] = r

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

    # Ensure gene calls carry an IMGT allele suffix (e.g. TRBV10-3 -> TRBV10-3*01).
    # Cell Ranger reports bare gene names, but tcrdist3's reference db is keyed by allele
    # and cannot map CDR1/CDR2 without it (empty clones -> crash). *01 is the safe default,
    # matching tcrdist3_matrix.py's own *00 -> *01 fallback.
    def add_allele(g):
        if g is None or pd.isna(g) or str(g) == '' or str(g) == 'None':
            return g
        s = str(g)
        return s if '*' in s else f"{s}*01"

    df['junction_aa'] = df['cdr3'].str.upper()
    df['v_call']      = (df['v_gene'] if 'v_gene' in df.columns else None)
    df['d_call']      = df['d_gene'] if 'd_gene' in df.columns else None
    df['j_call']      = df['j_gene'] if 'j_gene' in df.columns else None
    df['sequence']    = df['cdr3_nt'] if 'cdr3_nt' in df.columns else df['junction_aa']

    if 'v_gene' in df.columns:
        df['v_call'] = df['v_call'].apply(add_allele)
    if 'j_gene' in df.columns:
        df['j_call'] = df['j_call'].apply(add_allele)

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
        agg['sample']            = sample

        # Canonical schema + ORDER matching ANNOTATE_PROCESS so the shared engine's
        # positional sort/dedup ($1=junction_aa, $2=v_call) works. See sc_to_cdr3.py.
        out_cols = [
            'junction_aa', 'v_call', 'd_call', 'j_call',
            'duplicate_count', 'junction_aa_length', 'duplicate_frequency_percent',
            'sequence', 'sequence_id', 'junction', 'sample',
        ]
        out_path = f'bulk_samples/{sample}_bulk.tsv'
        agg[out_cols].to_csv(out_path, sep='\t', index=False)

        row = {'sample': sample, 'file': os.path.abspath(out_path)}
        smeta = meta_lookup.get(str(sample))
        for col in META_COLS:
            val = ''
            # Prefer the contigs, then the user sample sheet.
            if col in grp.columns:
                nn = grp[col].dropna()
                val = nn.iloc[0] if len(nn) > 0 else ''
            if (val == '' or pd.isna(val)) and smeta is not None and col in smeta.index:
                val = smeta[col]
            row[col] = '' if pd.isna(val) else val
        # Never leave patient_id empty (it drives PATIENT grouping) — fall back to sample.
        if row.get('patient_id', '') in ('', None) or pd.isna(row.get('patient_id', '')):
            row['patient_id'] = sample
        samplesheet_rows.append(row)

    ss = pd.DataFrame(samplesheet_rows)
    ss.to_csv('synthetic_samplesheet.csv', index=False)

    print(f"[Bridge 3] Created {len(samplesheet_rows)} bulk sample files in bulk_samples/")
    print(f"[Bridge 3] Synthetic samplesheet written to synthetic_samplesheet.csv")


if __name__ == '__main__':
    main()
