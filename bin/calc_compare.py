#!/usr/bin/env python3
"""
Description: this script calculates overlap measures between TCR repertoires

@author: Domenick Braccia
"""

import argparse
import pandas as pd
import numpy as np
from scipy.stats import entropy
import numpy as np
import csv

print('-- ENTERED calc_compare.py--')

# initialize parser
parser = argparse.ArgumentParser(description='Calculate clonality of a TCR repertoire')

# add arguments
parser.add_argument('-s', '--sample_table', 
                    metavar='sample_table', 
                    type=str, 
                    help='sample CSV file initially passed to nextflow run command')
parser.add_argument('-p', '--patient_table', 
                    metavar='patient_table', 
                    type=argparse.FileType('r'), 
                    help='patient CSV file initially passed to nextflow run command')

args = parser.parse_args() 

# Define functions
def jaccard_index(sample1, sample2):
    set1 = set(sample1)
    set2 = set(sample2)
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union

def sorensen_index(sample1, sample2):
    set1 = set(sample1)
    set2 = set(sample2)
    intersection = len(set1.intersection(set2))
    return 2 * intersection / (len(set1) + len(set2))

def morisita_horn_index(sample1, sample2):
    N1 = sum(sample1)
    N2 = sum(sample2)
    sum_n1i_n2i = sum([n1i * n2i for n1i, n2i in zip(sample1, sample2)])
    sum_n1i_sq = sum([n1i**2 for n1i in sample1])
    sum_n2i_sq = sum([n2i**2 for n2i in sample2])
    return 2 * sum_n1i_n2i / ((sum_n1i_sq + sum_n2i_sq) * (N1 + N2))


## Read in sample table CSV file
sample_table = pd.read_csv(args.sample_table, sep=',', header=0)
print('sample_table columns: \n')
print(sample_table.columns)

# Read in patient table CSV file
patient_table = pd.read_csv(args.patient_table, sep=',', header=0)
print('patient_table columns: \n')
print(patient_table.columns)

## Write dummy file named jaccard_amat.csv
with open('jaccard_amat.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['patient', 'sample', 'jaccard'])

## Write dummy file named sorensen_amat.csv
with open('sorensen_amat.csv', 'w') as f:
    writer = csv.writer(f)
    writer.writerow(['patient', 'sample', 'sorensen'])

