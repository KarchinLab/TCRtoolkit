#!/usr/bin/env python3

from cirro.helpers.preprocess_dataset import PreprocessDataset
import pandas as pd
import numpy as np

# 1. Get parameters from cirro pipeline call
ds = PreprocessDataset.from_running()
ds.logger.info("List of starting params")
ds.logger.info(ds.params)

ds.logger.info('checking ds.files')
files = ds.files
ds.logger.info(files.head())
ds.logger.info(files.columns)

# 2. Add samplesheet parameter and set equal to ds.samplesheet
ds.logger.info("Checking samplesheet parameter")
ds.logger.info(ds.samplesheet)
samplesheet = ds.samplesheet

# Replace local links with s3
samplesheet['file'] = files['file']

samplesheet.to_csv('samplesheet.csv', index=None)
ds.add_param("samplesheet", "samplesheet.csv")

ds.logger.info(ds.params)
