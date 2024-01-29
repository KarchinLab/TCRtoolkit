#!/usr/bin/env python3
"""
Description: this script calculates overlap measures between TCR repertoires

@author: Domenick Braccia
"""

import argparse
import pandas as pd
import numpy as np
import os
import csv
from scipy.stats import entropy

print('-- ENTERED calc_compare.py--')

# initialize parser
parser = argparse.ArgumentParser(description='Calculate clonality of a TCR repertoire')

# add arguments
parser.add_argument('-s', '--sample_utf8', 
                    metavar='sample_utf8', 
                    type=str, 
                    help='sample CSV file initially passed to nextflow run command')
parser.add_argument('-m', '--meta_data',
                    metavar='meta_data',
                    type=str,
                    help='metadata CSV file initially passed to nextflow run command')

args = parser.parse_args() 

## Read in sample table CSV file
## convert metadata to list
s = args.sample_utf8
sample_utf8 = pd.read_csv(args.sample_utf8, sep=',', header=0)
print('sample_utf8 looks like this: ' + str(sample_utf8))
print('sample_utf8 columns: \n')
print(sample_utf8.columns)

# Read in metadata table CSV file
meta_data = pd.read_csv(args.meta_data, sep=',', header=0)
print('meta_data looks like this: ' + str(meta_data))
print('meta_data columns: \n')
print(meta_data.columns)

# Import TCR count tables into dictionary of dataframes
file_paths = sample_utf8['file_path']
dfs = {}
for file_path in file_paths:
    # load data
    df = pd.read_csv(file_path, sep='\t', header=0)

    # Rename columns
    df = df.rename(columns={'count (templates/reads)': 'read_count', 'frequencyCount (%)': 'frequency'})
    dfs[file_path] = df

print('number of files in dfs: ' + str(len(dfs)))

## calculate the jaccard index between each sample pair in dfs and store in an nxn matrix and write to file
samples = dfs.keys()
for sample in samples:
    dfs[sample]

## calculate the sorensen index between each sample pair in dfs and store in an nxn matrix and write to file


## calculate the morisita index between each sample pair in dfs and store in an nxn matrix and write to file


## ========================================================================== ##

## Write dummy file named jaccard_amat.csv
with open('jaccard_amat.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['patient', 'sample', 'jaccard'])

## Write dummy file named sorensen_amat.csv
with open('sorensen_amat.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['patient', 'sample', 'sorensen'])

## Write dummy file named morisita_amat.csv
with open('morisita_amat.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['patient', 'sample', 'morisita'])
