#!/usr/bin/env python3
"""
Bridge: SC_TO_CDR3
Converts TCELL_INTEGRATION per-cell export_cells.tsv directly into
CDR3b/TRBV/TRBJ/counts/sample format, bypassing ANNOTATE_PROCESS.

Cells are grouped per (sample, CDR3b, TRBV, TRBJ); the cell count
becomes the clonotype 'counts' field.

Outputs:
    per_sample/{sample}_cdr3.tsv  — one file per sample
    concat_cdr3.tsv               — all samples concatenated (input for ANNOTATE_SORT_CDR3)
"""

import sys
import os
import pandas as pd


def main():
    if len(sys.argv) < 2:
        print("Usage: sc_to_cdr3.py <export_cells.tsv>")
        sys.exit(1)

    export_cells_file = sys.argv[1]
    df = pd.read_csv(export_cells_file, sep='\t', low_memory=False)

    # Resolve CDR3b
    if 'cdr3b' in df.columns:
        df['CDR3b'] = df['cdr3b']
    elif 'junction_aa' in df.columns:
        df['CDR3b'] = df['junction_aa']
    else:
        raise ValueError("Export file must contain 'cdr3b' or 'junction_aa'.")

    # Resolve TRBV / TRBJ
    df['TRBV'] = df['trbv'] if 'trbv' in df.columns else None
    df['TRBJ'] = df['trbj'] if 'trbj' in df.columns else None

    if 'sample' not in df.columns:
        raise ValueError("Export file must contain a 'sample' column.")

    # Keep only cells with a valid beta CDR3
    df = df[df['CDR3b'].notna() & (df['CDR3b'] != '') & (df['CDR3b'] != 'None')]
    print(f"[SC_TO_CDR3] {len(df)} cells with valid CDR3b across {df['sample'].nunique()} sample(s).")

    os.makedirs('per_sample', exist_ok=True)
    all_frames = []

    for sample, grp in df.groupby('sample'):
        agg = (
            grp.groupby(['CDR3b', 'TRBV', 'TRBJ'], dropna=False)
               .size()
               .reset_index(name='counts')
        )
        agg['sample'] = sample
        agg = agg[['CDR3b', 'TRBV', 'TRBJ', 'counts', 'sample']]

        out_path = f'per_sample/{sample}_cdr3.tsv'
        agg.to_csv(out_path, sep='\t', index=False)
        all_frames.append(agg)
        print(f"[SC_TO_CDR3]   {sample}: {len(agg)} clonotypes.")

    if not all_frames:
        raise ValueError("[SC_TO_CDR3] No clonotypes found after filtering.")

    concat = pd.concat(all_frames, ignore_index=True)
    concat.to_csv('concat_cdr3.tsv', sep='\t', index=False)
    print(f"[SC_TO_CDR3] concat_cdr3.tsv written ({len(concat)} rows total).")


if __name__ == '__main__':
    main()
