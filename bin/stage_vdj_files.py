#!/usr/bin/env python3
"""
Stages Cirro/S3-resolved VDJ files (contigs/clonotypes/metrics, staged
individually by Nextflow into numbered subdirectories - Nextflow can't stage
a whole directory straight from S3) back into a reconstructed local outs/
directory per sample, and rewrites the sample_sheet's path column to point at
it.

Positional convention: subdirectory N under each --*-dir corresponds to row N
(0-indexed) of --sample-sheet, matching the order .cirro_singlecell_*/
preprocess.py emitted the per-sample file lists in. A subdirectory containing
only a file named NO_FILE means nothing was staged for that sample at that
position (local/on-prem runs, or no match found upstream) - that row's `path`
column is left exactly as-is.

Usage:
    stage_vdj_files.py --sample-sheet sample_sheet.csv \\
        --contigs-dir staged_vdj/contigs_input \\
        --clonotypes-dir staged_vdj/clonotypes_input \\
        --metrics-dir staged_vdj/metrics_input \\
        --out-sample-sheet resolved_sample_sheet.csv
"""

import argparse
import csv
import glob
import os
import shutil

NO_FILE_MARKER = "NO_FILE"


def staged_file_for(stage_dir_prefix, index):
    """
    Return the path Nextflow staged into <stage_dir_prefix><index+1>/ (the
    'contigs_input*/*' stageAs pattern numbers subdirectories from 1), or
    None if that position has no real file - either the NO_FILE sentinel, or
    the subdirectory doesn't exist because fewer items were staged than rows.
    """
    candidates = glob.glob(f"{stage_dir_prefix}{index + 1}/*")
    if not candidates:
        return None
    fp = candidates[0]
    return None if os.path.basename(fp) == NO_FILE_MARKER else fp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-sheet", required=True)
    ap.add_argument("--contigs-dir", required=True)
    ap.add_argument("--clonotypes-dir", required=True)
    ap.add_argument("--metrics-dir", required=True)
    ap.add_argument("--out-sample-sheet", required=True)
    args = ap.parse_args()

    with open(args.sample_sheet, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)

    for i, row in enumerate(rows):
        contigs_fp = staged_file_for(args.contigs_dir, i)
        clonotypes_fp = staged_file_for(args.clonotypes_dir, i)
        metrics_fp = staged_file_for(args.metrics_dir, i)

        if contigs_fp is None and clonotypes_fp is None and metrics_fp is None:
            # Nothing staged for this sample - local/on-prem run, or no
            # match found upstream. Leave the row's own path untouched.
            continue

        outs_dir = os.path.join("resolved_outs", f"sample_{i}", "outs")
        os.makedirs(outs_dir, exist_ok=True)
        for fp, real_name in (
            (contigs_fp, "filtered_contig_annotations.csv"),
            (clonotypes_fp, "clonotypes.csv"),
            (metrics_fp, "metrics_summary.csv"),
        ):
            if fp is not None:
                shutil.copy(fp, os.path.join(outs_dir, real_name))

        print(f"[stage_vdj_files] sample '{row['sample']}': staged VDJ files into {outs_dir}")
        row["path"] = os.path.join("resolved_outs", f"sample_{i}")

    with open(args.out_sample_sheet, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
