#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
import numpy as np
from pathlib import Path

# TCRtoolkit single-cell mode needs a 'sample_sheet.csv' with a 'path' column
# pointing at each sample's Cell Ranger VDJ outs/ DIRECTORY, not a single
# file - see modules/scratch/VDJ_QC/VDJ_QC_analysis.qmd's column contract
# (sample_sheet_sample_col="sample", sample_sheet_path_col="path"). Adapted
# from break-through-cancer/staple's .cirro/preprocess.py, which solves the
# same file-vs-directory mismatch for its own per-sample directory inputs.
SAMPLESHEET_REQUIRED_COLUMNS = ("sample", "path")


def samplesheet_from_files(ds):
    files = ds.files
    ds.logger.info(f"found files in ds.files: {files.to_dict()}")

    # Cirro's per-file listing gives one row per file under each sample's
    # directory; the sample's directory is that file's parent. Path()
    # collapses "s3://" to "s3:/", so restore the double slash afterward.
    files['path'] = files['file'].apply(lambda x: str(Path(x).parent).replace('s3:/', 's3://'))
    files = files[['sample', 'path']].drop_duplicates()

    return pd.merge(ds.samplesheet, files, on='sample', how='left')


def samplesheet_from_params(ds):
    # Fallback for datasets registered as one directory-shaped input per
    # sample (no flat per-file listing in ds.files).
    return pd.DataFrame({
        'sample': [x['name'] for x in ds.metadata['inputs']],
        'path': [x['dataPath'] for x in ds.metadata['inputs']],
    })


def prepare_sample_sheet(ds):
    sample_sheet = samplesheet_from_files(ds)

    if sample_sheet.empty:
        ds.logger.warning("No files found in dataset. Preparing sample_sheet from params.")
        sample_sheet = samplesheet_from_params(ds)
        if sample_sheet.empty:
            raise ValueError("No files found in dataset and unable to prepare sample_sheet from params.")

    for colname in SAMPLESHEET_REQUIRED_COLUMNS:
        if colname not in sample_sheet.columns:
            ds.logger.warning(f"sample_sheet is missing required column '{colname}'. Populating with NaN.")
            sample_sheet[colname] = np.nan

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
