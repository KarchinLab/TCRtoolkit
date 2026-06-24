#!/usr/bin/env python3
"""
Bridge 2: CLUSTER_TO_SC  (optional)
Maps TCRtoolkit bulk cluster assignments (CDR3b-level) onto single cells
using CDR3b as the join key, producing per-cell TSVs compatible with
SCRATCH-TCR's Consensus Clustering module inputs.

Usage:
    cluster_to_sc.py <export_cells.tsv> <giana_clusters> <gliph_clusters>

Pass 'NO_FILE' for any cluster file that is not available.

Outputs (written only when input is not NO_FILE):
    giana_export_cells.tsv   — per-cell with giana_cluster column
    gliph2_export_cells.tsv  — per-cell with gliph2_cluster column
"""

import sys
import pandas as pd


NO_FILE_SENTINEL = 'NO_FILE'


def extract_beta_cdr3(ctaa):
    """Extract beta chain CDR3 from CTaa field (format: 'alpha;beta')."""
    if pd.isna(ctaa):
        return None
    parts = str(ctaa).split(';')
    return parts[1].strip().upper() if len(parts) > 1 else None


def map_clusters(cells: pd.DataFrame, cluster_file: str, cluster_col_name: str,
                 out_file: str) -> None:
    """
    Join CDR3b-level cluster assignments from a TCRtoolkit output file onto
    per-cell data and write the result.

    The cluster file is expected to have:
      - First column: CDR3b amino acid sequence (may have a comment header)
      - 'cluster' column (or second column if 'cluster' is absent)
    """
    try:
        cl = pd.read_csv(cluster_file, sep='\t', comment='#')
    except Exception as e:
        print(f"[Bridge 2] Warning: could not read {cluster_file}: {e}")
        return

    cdr3_col    = cl.columns[0]
    cluster_col = 'cluster' if 'cluster' in cl.columns else cl.columns[1]

    cl_map = (
        cl[[cdr3_col, cluster_col]]
        .rename(columns={cdr3_col: 'CDR3b', cluster_col: cluster_col_name})
        .drop_duplicates('CDR3b')
    )
    cl_map['CDR3b'] = cl_map['CDR3b'].str.upper()

    merged = cells.merge(cl_map, on='CDR3b', how='left')
    merged.to_csv(out_file, sep='\t', index=False)
    assigned = merged[cluster_col_name].notna().sum()
    print(f"[Bridge 2] {out_file}: {assigned}/{len(merged)} cells assigned a {cluster_col_name}")


def main():
    if len(sys.argv) < 4:
        print("Usage: cluster_to_sc.py <export_cells.tsv> <giana_clusters> <gliph_clusters>")
        sys.exit(1)

    export_cells_file = sys.argv[1]
    giana_file        = sys.argv[2]
    gliph_file        = sys.argv[3]

    cells = pd.read_csv(export_cells_file, sep='\t', low_memory=False)

    if 'CTaa' in cells.columns:
        cells['CDR3b'] = cells['CTaa'].apply(extract_beta_cdr3)
    elif 'junction_aa' in cells.columns:
        cells['CDR3b'] = cells['junction_aa'].str.upper()
    else:
        raise ValueError("export_cells.tsv must contain 'CTaa' or 'junction_aa'")

    if giana_file and giana_file != NO_FILE_SENTINEL:
        map_clusters(cells, giana_file, 'giana_cluster', 'giana_export_cells.tsv')

    if gliph_file and gliph_file != NO_FILE_SENTINEL:
        map_clusters(cells, gliph_file, 'gliph2_cluster', 'gliph2_export_cells.tsv')


if __name__ == '__main__':
    main()
