#!/usr/bin/env python3
"""
Bridge 1: SC_TO_BULK
Converts SCRATCH-TCR per-cell export_cells.tsv → per-sample AIRR-format TSVs.

Usage:
    sc_to_bulk.py <export_cells.tsv> <sample_col> [<meta_cols_csv>]

Outputs:
    bulk_samples/<sample>_bulk.tsv   — one AIRR TSV per sample
    synthetic_samplesheet.csv        — TCRtoolkit-compatible samplesheet

Column resolution priority for junction_aa / v_call / j_call:
  1. Dedicated columns: cdr3b, trbv, trbj  (TCELL_INTEGRATION export format)
  2. CTaa / CTgene parsing                  (legacy scRepertoire semicolon format)
  3. Pre-existing junction_aa column        (already AIRR format)
"""

import sys
import os
import pandas as pd


def extract_beta_ctaa(field):
    """Extract beta CDR3 from CTaa.
    Handles both:
      - semicolon-separated: 'ALPHA;BETA'
      - pipe + label format: 'A:ALPHA|B:BETA'
    """
    if pd.isna(field) or not field:
        return None
    s = str(field)
    # pipe + label format (A:...|B:...)
    if '|' in s:
        for part in s.split('|'):
            part = part.strip()
            if part.upper().startswith('B:'):
                return part[2:] or None
        return None
    # semicolon-separated format
    parts = s.split(';')
    return parts[1] if len(parts) > 1 else None


def extract_beta_ctgene(field):
    """Extract beta gene string from CTgene (same split logic as CTaa)."""
    return extract_beta_ctaa(field)


def extract_trbv(ctgene_beta):
    if not ctgene_beta:
        return None
    genes = ctgene_beta.split('.')
    return genes[0] if genes else None


def extract_trbj(ctgene_beta):
    if not ctgene_beta:
        return None
    genes = ctgene_beta.split('.')
    return genes[-1] if len(genes) > 1 else None


def main():
    if len(sys.argv) < 3:
        print("Usage: sc_to_bulk.py <export_cells.tsv> <sample_col> [meta_cols_csv]")
        sys.exit(1)

    export_cells_file = sys.argv[1]
    sample_col        = sys.argv[2]
    meta_cols         = [c for c in sys.argv[3].split(',') if c] if len(sys.argv) > 3 else []

    df = pd.read_csv(export_cells_file, sep='\t', low_memory=False)

    # ── Resolve junction_aa (beta CDR3) ──────────────────────────────────────
    # Priority: dedicated cdr3b col > CTaa parsing > existing junction_aa col
    if 'cdr3b' in df.columns:
        df['junction_aa'] = df['cdr3b'].where(df['cdr3b'].notna() & (df['cdr3b'] != ''))
        print("[Bridge 1] Using 'cdr3b' column for junction_aa.")
    elif 'CTaa' in df.columns:
        df['junction_aa'] = df['CTaa'].apply(extract_beta_ctaa)
        print("[Bridge 1] Parsed junction_aa from CTaa column.")
    elif 'junction_aa' in df.columns:
        print("[Bridge 1] Using pre-existing junction_aa column.")
    else:
        raise ValueError("export_cells.tsv must contain 'cdr3b', 'CTaa', or 'junction_aa'.")

    # ── Resolve v_call / j_call ───────────────────────────────────────────────
    if 'trbv' in df.columns:
        df['v_call'] = df['trbv'].where(df['trbv'].notna() & (df['trbv'] != ''))
        print("[Bridge 1] Using 'trbv' column for v_call.")
    elif 'CTgene' in df.columns:
        beta_gene    = df['CTgene'].apply(extract_beta_ctgene)
        df['v_call'] = beta_gene.apply(extract_trbv)
    else:
        df['v_call'] = None

    if 'trbj' in df.columns:
        df['j_call'] = df['trbj'].where(df['trbj'].notna() & (df['trbj'] != ''))
        print("[Bridge 1] Using 'trbj' column for j_call.")
    elif 'CTgene' in df.columns:
        beta_gene    = df['CTgene'].apply(extract_beta_ctgene)
        df['j_call'] = beta_gene.apply(extract_trbj)
    else:
        df['j_call'] = None

    # ── Keep only cells with a valid beta CDR3 ────────────────────────────────
    df = df[df['junction_aa'].notna() & (df['junction_aa'] != 'None') & (df['junction_aa'] != '')]
    print(f"[Bridge 1] {len(df)} cells retained after filtering for valid junction_aa.")

    # ── Resolve sample column ─────────────────────────────────────────────────
    if sample_col not in df.columns:
        if 'sample' in df.columns:
            print(f"[Bridge 1] WARNING: sample_col '{sample_col}' not found; falling back to 'sample'.")
            sample_col = 'sample'
        else:
            raise ValueError(f"sample_col '{sample_col}' not found in export_cells.tsv. "
                             f"Available columns: {list(df.columns)}")

    os.makedirs('bulk_samples', exist_ok=True)
    samplesheet_rows = []

    for sample, grp in df.groupby(sample_col):
        agg = (
            grp.groupby(['junction_aa', 'v_call', 'j_call'])
               .size()
               .reset_index(name='duplicate_count')
        )
        total = agg['duplicate_count'].sum()
        agg['duplicate_frequency_percent'] = (agg['duplicate_count'] / total * 100).round(6)
        agg['sequence_id'] = agg['junction_aa']
        agg['sequence']    = agg['junction_aa']  # placeholder; NT not available from per-cell table

        out_cols = [
            'sequence_id', 'junction_aa', 'v_call', 'j_call',
            'duplicate_count', 'duplicate_frequency_percent', 'sequence'
        ]
        out_path = f'bulk_samples/{sample}_bulk.tsv'
        agg[out_cols].to_csv(out_path, sep='\t', index=False)

        row = {'sample': sample, 'file': out_path}
        for col in meta_cols:
            if col in grp.columns:
                row[col] = grp[col].iloc[0]
        samplesheet_rows.append(row)

    ss = pd.DataFrame(samplesheet_rows)
    ss.to_csv('synthetic_samplesheet.csv', index=False)

    print(f"[Bridge 1] Created {len(samplesheet_rows)} bulk sample files in bulk_samples/")
    print(f"[Bridge 1] Synthetic samplesheet written to synthetic_samplesheet.csv")


if __name__ == '__main__':
    main()
