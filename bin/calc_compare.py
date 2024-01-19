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

