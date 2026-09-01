#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
import numpy as np
from pathlib import Path

# TCRtoolkit single-cell mode needs a 'sample_sheet.csv' with a 'path' column
# pointing at each sample's Cell Ranger VDJ outs/ DIRECTORY - see
# modules/scratch/VDJ_QC/VDJ_QC_analysis.qmd's column contract
# (sample_sheet_sample_col="sample", sample_sheet_path_col="path").
#
# Cirro's dataset-file abstraction doesn't stage whole directories, only
# individual files (confirmed against WangLab-ComputationalBiology/SCRATCH-QC's
# own .cirro/qc/preprocess.py, a sibling pipeline solving the same problem: it
# matches specific files by name via ds.files rather than referencing a
# directory). So this profile finds the 3 specific files VDJ_QC actually reads
# out of each sample's outs/ directory (see VDJ_QC_analysis.qmd's
# detect_vdj_file()) and exposes them as separate, sample-ordered pipeline
# params. modules/scratch/VDJ_QC/main.nf stages them and reconstructs a local
# outs/ directory per sample, rewriting sample_sheet.csv's path column to
# point at the staged copy before VDJ_QC's R code (unchanged) reads it.
SAMPLESHEET_REQUIRED_COLUMNS = ("sample", "path")

# Maps each file VDJ_QC's detect_vdj_file() looks for to the Nextflow param
# modules/scratch/VDJ_QC/main.nf reads it back from.
VDJ_FILE_PARAMS = {
    "filtered_contig_annotations.csv": "contigs_files",
    "clonotypes.csv": "clonotypes_files",
    "metrics_summary.csv": "metrics_files",
}


def samplesheet_from_files(ds):
    files = ds.files
    ds.logger.info(f"found files in ds.files: {files.to_dict()}")

    # Cirro's per-file listing gives one row per file under each sample's
    # directory; the sample's directory is that file's parent. Path()
    # collapses "s3://" to "s3:/", so restore the double slash afterward.
    files = files.copy()
    files['path'] = files['file'].apply(lambda x: str(Path(x).parent).replace('s3:/', 's3://'))
    files = files[['sample', 'path']].drop_duplicates()

    return pd.merge(ds.samplesheet, files, on='sample', how='left')


def samplesheet_from_params(ds):
    # Fallback for datasets registered as one directory-shaped input per
    # sample (no flat per-file listing in ds.files). Not paired with
    # prepare_vdj_file_params() - a genuine directory-shaped registration
    # gives one real, already-resolvable path per sample directly.
    return pd.DataFrame({
        'sample': [x['name'] for x in ds.metadata['inputs']],
        'path': [x['dataPath'] for x in ds.metadata['inputs']],
    })


def prepare_vdj_file_params(ds, samples):
    """
    Match each sample to its filtered_contig_annotations.csv / clonotypes.csv /
    metrics_summary.csv via ds.files's flat per-file listing. Sets one
    Nextflow param per filename (VDJ_FILE_PARAMS), each a comma-separated,
    sample-ordered list of S3 paths - empty string where a sample has no
    match, so list position always lines up 1:1 with `samples` even when a
    file is missing for one sample. Uses split(',', -1) semantics on the
    Nextflow side, so trailing empty entries can't be silently dropped.
    """
    files = ds.files

    for fname, param_name in VDJ_FILE_PARAMS.items():
        hits = files[files['file'].str.endswith(fname)]
        per_sample = []
        for sample in samples:
            sample_hits = sorted(hits.loc[hits['sample'] == sample, 'file'].tolist())
            if not sample_hits:
                ds.logger.warning(
                    f"No {fname} found for sample '{sample}' via ds.files - "
                    "VDJ_QC will fall back to sample_sheet's own path column for this sample."
                )
                per_sample.append('')
            else:
                if len(sample_hits) > 1:
                    ds.logger.warning(f"Multiple {fname} matches for sample '{sample}', using the first: {sample_hits}")
                per_sample.append(sample_hits[0])
        ds.add_param(param_name, ','.join(per_sample))


def prepare_sample_sheet(ds):
    sample_sheet = samplesheet_from_files(ds)

    # A left merge against an empty ds.files still returns every ds.samplesheet
    # row - just with 'path' as NaN - so sample_sheet.empty alone can't detect
    # "no per-file listing to merge against". Check the column directly.
    if sample_sheet.empty or sample_sheet['path'].isna().all():
        ds.logger.warning("No files found in dataset. Preparing sample_sheet from params.")
        sample_sheet = samplesheet_from_params(ds)
        if sample_sheet.empty:
            raise ValueError("No files found in dataset and unable to prepare sample_sheet from params.")
    else:
        prepare_vdj_file_params(ds, sample_sheet['sample'].tolist())

    for colname in SAMPLESHEET_REQUIRED_COLUMNS:
        if colname not in sample_sheet.columns:
            ds.logger.warning(f"sample_sheet is missing required column '{colname}'. Populating with NaN.")
            sample_sheet[colname] = np.nan

    try:
        # The PATIENT step (GIANA + GLIPH2) pools samples by params.patient_col. If that column
        # is absent every sample becomes its own patient and the clustering is silently wrong -
        # no error, just useless results - so say so loudly here. samplesheet_from_params()'s
        # fallback frame in particular carries only sample/path.
        # ds.params may be a plain dict or a params object depending on cirro version, and a
        # warning must never be the thing that breaks preprocessing - so degrade to the default.
        try:
            patient_col = dict(ds.params).get("patient_col") or "patient_id"
        except Exception:
            patient_col = "patient_id"
        if patient_col not in sample_sheet.columns:
            ds.logger.warning(
                f"sample_sheet has no '{patient_col}' column. GIANA and GLIPH2 pool samples per "
                "patient, so each sample will be treated as its own patient and cross-sample "
                "clustering within a patient will be lost. Add the column to the dataset "
                "samplesheet, or set 'Patient column' to one that exists."
            )
        elif sample_sheet[patient_col].isna().any():
            missing = sample_sheet.loc[sample_sheet[patient_col].isna(), "sample"].tolist()
            ds.logger.warning(f"Samples with no {patient_col}: {missing}. These will not pool with any patient.")
        else:
            n_pat = sample_sheet[patient_col].nunique()
            ds.logger.info(f"{len(sample_sheet)} samples across {n_pat} patients (by '{patient_col}')")
    except Exception as e:
        ds.logger.warning(f"patient-column check skipped: {e}")

    sample_sheet.to_csv('sample_sheet.csv', index=None)
    ds.add_param('sample_sheet', 'sample_sheet.csv')
    ds.logger.info(sample_sheet.to_dict())


def main():
    ds = PreprocessDataset.from_running()
    ds.logger.info("List of starting params")
    ds.logger.info(ds.params)

    prepare_sample_sheet(ds)

    ds.logger.info(ds.params)


if __name__ == "__main__":
    main()
