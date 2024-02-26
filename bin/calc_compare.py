#!/usr/bin/env python3
"""
Description: this script calculates overlap measures between TCR repertoires

@author: Domenick Braccia
"""

import argparse
import pandas as pd
import numpy as np
import os
import sys
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
parser.add_argument('-p', '--project_dir',
                    metavar='project_dir',
                    type=str,
                    help='path to project directory')

args = parser.parse_args() 

## Import project directory path
project_dir = args.project_dir
sys.path.append(project_dir + '/bin/')
from utils import TicTocGenerator, tic, toc, jaccard_index, sorensen_index, morisita_horn_index
TicToc = TicTocGenerator()

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
samples = list(dfs.keys())
jaccard_mat = np.zeros((len(samples), len(samples)))
for i, sample1 in enumerate(samples):
    for j, sample2 in enumerate(samples):
        # calculate jaccard index
        value = jaccard_index(dfs[sample1]['aminoAcid'], dfs[sample2]['aminoAcid'])
        # store in numpy array
        jaccard_mat[i, j] = value

# define column and index names
sample_names= [os.path.basename(sample).split('.')[0] for sample in samples]
jaccard_df = pd.DataFrame(jaccard_mat, columns=sample_names, index=sample_names)

# save jacard_df to csv
jaccard_df.to_csv('jaccard_mat.csv', index=True, header=True)

## calculate the sorensen index between each sample pair in dfs and store in an nxn matrix and write to file
sorensen_mat = np.zeros((len(samples), len(samples)))
for i, sample1 in enumerate(samples):
    for j, sample2 in enumerate(samples):
        # calculate sorensen index
        value = sorensen_index(dfs[sample1]['aminoAcid'], dfs[sample2]['aminoAcid'])
        # store in numpy array
        sorensen_mat[i, j] = value

# define column and index names
sorensen_df = pd.DataFrame(sorensen_mat, columns=sample_names, index=sample_names)

# save sorensen_df to csv
sorensen_df.to_csv('sorensen_mat.csv', index=True, header=True)

## calculate the morisita index between each sample pair in dfs and store in an nxn matrix and write to file
morisita_mat = np.zeros((len(samples), len(samples)))
for i in range(len(samples)):
    for j in range(i+1):
        print('- on samples: [' + str(i) + ', ' + str(j) + ']')
        # calculate morisita index
        value = morisita_horn_index(dfs, samples[i], samples[j])
        # store in numpy array
        morisita_mat[i, j] = value

# Copy the lower triangle to the upper triangle
morisita_mat = morisita_mat + morisita_mat.T - np.diag(morisita_mat.diagonal())

# morisita_mat = np.zeros((len(samples), len(samples)))
# for i, sample1 in enumerate(samples):
#     for j, sample2 in enumerate(samples):
#         print('on samples: [' + str(i) + ', ' + str(j) + ']')
#         # calculate morisita index
#         value = morisita_horn_index(dfs, sample1, sample2)
#         # store in numpy array
#         morisita_mat[i, j] = value

# define column and index names
morisita_df = pd.DataFrame(morisita_mat, columns=sample_names, index=sample_names)

# save morisita_df to csv
morisita_df.to_csv('morisita_mat.csv', index=True, header=True)

## ========================================================================== ##
