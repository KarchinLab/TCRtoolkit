#!/usr/bin/env python3
"""
Pseudobulk QC gate — per-sample metrics + gene-family usage.

Computes a lightweight QC summary for a single pseudobulk-derived TCR sample
table produced by either single-cell → bulk bridge:

  - SC_TO_BULK  (AIRR-ish):  junction_aa, v_call, j_call, duplicate_count
  - SC_TO_CDR3            :  CDR3b,       TRBV,   TRBJ,   counts

Metrics (intentionally minimal — diversity metrics such as Shannon/Simpson are
NOT computed here because the REPERTOIRE module already reports them):

  n_clones : number of unique clonotype rows
  n_cells  : total cells (sum of the per-clone count column)

A sample PASSes the gate when n_clones >= min_clones AND n_cells >= min_cells.

Also emits non-redundant gene-family usage (V-family clone counts, TRBV1..30)
for reporting.

Usage:
    pseudobulk_qc.py <sample_tsv> <sample_name> <min_clones> <min_cells>

Outputs (cwd):
    qc_<sample>.csv        sample,n_clones,n_cells,min_clones,min_cells,pass
    vfamily_<sample>.csv   sample,TRBV1,...,TRBV30
"""

import sys
import re
import pandas as pd


def resolve_col(df, *candidates):
    """Return the first candidate column present in df, else None."""
    for c in candidates:
        if c in df.columns:
            return c
    return None


def extract_trbv_family(allele):
    """Map a V allele (e.g. 'TRBV20-1*01') to its family ('TRBV20')."""
    if pd.isna(allele):
        return None
    m = re.match(r'(TRBV)(\d+)', str(allele))
    return f"{m.group(1)}{m.group(2)}" if m else None


def main():
    if len(sys.argv) < 5:
        sys.exit("Usage: pseudobulk_qc.py <sample_tsv> <sample_name> <min_clones> <min_cells>")

    sample_tsv = sys.argv[1]
    sample     = sys.argv[2]
    min_clones = int(sys.argv[3])
    min_cells  = int(sys.argv[4])

    df = pd.read_csv(sample_tsv, sep='\t', low_memory=False)

    # Auto-detect the count column (AIRR 'duplicate_count' vs SC_TO_CDR3 'counts')
    count_col = resolve_col(df, 'duplicate_count', 'counts', 'count')
    v_col     = resolve_col(df, 'v_call', 'TRBV', 'trbv')

    n_clones = int(len(df))
    if count_col is not None:
        n_cells = int(pd.to_numeric(df[count_col], errors='coerce').fillna(0).sum())
    else:
        # No count column → treat each row as a single cell.
        n_cells = n_clones

    passed = (n_clones >= min_clones) and (n_cells >= min_cells)

    pd.DataFrame([{
        'sample': sample,
        'n_clones': n_clones,
        'n_cells': n_cells,
        'min_clones': min_clones,
        'min_cells': min_cells,
        'pass': 'PASS' if passed else 'FAIL',
    }]).to_csv(f'qc_{sample}.csv', index=False)

    # ── V gene family usage (clone counts per TRBV family) ────────────────────
    all_fams = [f'TRBV{i}' for i in range(1, 31)]
    if v_col is not None:
        fam_counts = df[v_col].apply(extract_trbv_family).value_counts(dropna=True)
    else:
        fam_counts = pd.Series(dtype=int)

    row = {'sample': sample}
    for fam in all_fams:
        row[fam] = int(fam_counts.get(fam, 0))
    pd.DataFrame([row]).to_csv(f'vfamily_{sample}.csv', index=False)

    # Human-readable log line
    print(f"[Pseudobulk QC] {sample}: n_clones={n_clones}, n_cells={n_cells} -> "
          f"{'PASS' if passed else 'FAIL'} (min_clones={min_clones}, min_cells={min_cells})")


if __name__ == '__main__':
    main()
