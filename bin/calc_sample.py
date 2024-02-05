#!/usr/bin/env python3
"""
Description: this script calculates the clonality of a TCR repertoire

@author: Domenick Braccia
@contributor: elhanaty
"""

## import packages
import argparse
import pandas as pd
import numpy as np
from scipy.stats import entropy
import numpy as np
import csv

# print('-- ENTERED calc_clonality.py--')

# initialize parser
parser = argparse.ArgumentParser(description='Calculate clonality of a TCR repertoire')

# add arguments
parser.add_argument('-s', '--sample_meta', 
                    metavar='sample_meta', 
                    type=str, 
                    help='sample metadata passed in through samples CSV file')
parser.add_argument('-c', '--count_table', 
                    metavar='count_table', 
                    type=argparse.FileType('r'), 
                    help='counts file in TSV format')
parser.add_argument('-m', '--meta_data',
                    metavar='meta_data',
                    type=argparse.FileType('r'),
                    help='metadata file in CSV format')

args = parser.parse_args() 

## convert metadata to list
s = args.sample_meta
sample_meta = args.sample_meta[1:-1].split(', ')
# print('sample_meta looks like this: ' + str(sample_meta))

# Read in the counts file
counts = pd.read_csv(args.count_table, sep='\t', header=0)
counts = counts.rename(columns={'count (templates/reads)': 'read_count', 'frequencyCount (%)': 'frequency'})
# print('counts columns: \n')
# print(counts.columns)

def calc_sample_stats(sample_meta, counts):
    """Calculate sample level statistics of TCR repertoire."""

    ## first pass stats
    clone_counts = counts['read_count']
    clone_entropy = entropy(clone_counts, base=2)
    num_clones = len(clone_counts)
    num_TCRs = sum(clone_counts)
    clonality = 1 - clone_entropy / np.log2(num_clones)
    simpson_index = sum(clone_counts**2)/(num_TCRs**2)
    simpson_index_corrected = sum(clone_counts*(clone_counts-1))/(num_TCRs*(num_TCRs-1))

    ## tcr productivity stats
    clone_prod = counts['sequenceStatus']
    # print('clone_prod looks like this: ' + str(clone_prod))

    # count number of productive clones
    num_in = sum(clone_prod == 'In')
    num_out = sum(clone_prod == 'Out')
    num_stop = sum(clone_prod == 'Stop')
    pct_prod = num_in / num_clones
    pct_out = num_out / num_clones
    pct_stop = num_stop / num_clones
    pct_nonprod = pct_out + pct_stop

    ## cdr3 info
    cdr3_lens = counts['cdr3Length']
    cdr3_avg_len = np.mean(cdr3_lens)

    ## Calculate convergence for each T cell receptor
    aas = counts[counts.aminoAcid.notnull()].aminoAcid.unique()
    dict_df = {}
    for aa in aas:
        dict_df[aa] = {'counts': counts[counts.aminoAcid == aa]}
        # append key value pair to dict_df[aa] with key convergence equal to the number of rows in counts
        dict_df[aa]['convergence'] = len(counts[counts.aminoAcid == aa])

    ## calculate the number of covergent TCRs for each sample
    num_convergent = 0
    for aa in aas:
      if dict_df[aa]['convergence'] > 1:
        num_convergent += 1    
    
    ## calculate ratio of convergent TCRs to total TCRs
    ratio_convergent = num_convergent/len(aas)

    ## add in patient meta_data such as responder status to sample_stats.csv
    # read in metadata file
    meta_data = pd.read_csv(args.meta_data, sep=',', header=0)

    # filter out metadata for the current sample
    current_meta = meta_data[meta_data['patient_id'] == sample_meta[1]]

    # write above values to csv file
    with open('sample_stats.csv', 'w') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([sample_meta[0], sample_meta[1], sample_meta[2], sample_meta[3],
                         str(current_meta.treatment.iloc[0]), str(current_meta.response.iloc[0]), str(current_meta.clinical_data.iloc[0]),
                         num_clones, num_TCRs, simpson_index, simpson_index_corrected, clonality,
                         num_in, num_out, num_stop, pct_prod, pct_out, pct_stop, pct_nonprod,
                         cdr3_avg_len, num_convergent, ratio_convergent])
        
    # store v_family gene usage in a dataframe
    v_family = counts['vFamilyName'].value_counts(dropna=False).to_frame().T.sort_index(axis=1)
    d_family = counts['dFamilyName'].value_counts(dropna=False).to_frame().T.sort_index(axis=1)
    j_family = counts['jFamilyName'].value_counts(dropna=False).to_frame().T.sort_index(axis=1)

    # generate a list of all possible columns names from TCRBV01-TCRBV30
    all_v_fam = ['TCRBV{:02d}'.format(i) for i in range(1, 31)]

    # generate a list of all possible columns names from TCRBD01-TCRBD02
    all_d_fam = ['TCRBD{:02d}'.format(i) for i in range(1, 3)]

    # generate a list of all possible columns names from TCRBJ01-TCRBJ02
    all_j_fam = ['TCRBJ{:02d}'.format(i) for i in range(1, 3)]

    # add missing columns to v_family dataframe by reindexing
    v_family_reindex = v_family.reindex(columns=all_v_fam, fill_value=0)
    d_family_reindex = d_family.reindex(columns=all_d_fam, fill_value=0)
    j_family_reindex = j_family.reindex(columns=all_j_fam, fill_value=0)

    # add sample_meta columns to v_family_reindex and make them the first three columns
    v_family_reindex.insert(0, 'origin', sample_meta[3])
    v_family_reindex.insert(0, 'timepoint', sample_meta[2])
    v_family_reindex.insert(0, 'patient_id', sample_meta[1])
    d_family_reindex.insert(0, 'origin', sample_meta[3])
    d_family_reindex.insert(0, 'timepoint', sample_meta[2])
    d_family_reindex.insert(0, 'patient_id', sample_meta[1])
    j_family_reindex.insert(0, 'origin', sample_meta[3])
    j_family_reindex.insert(0, 'timepoint', sample_meta[2])
    j_family_reindex.insert(0, 'patient_id', sample_meta[1])

    # Write v_family_reindex to csv file with no header and no index
    v_family_reindex.to_csv('v_family.csv', header=False, index=False)
    d_family_reindex.to_csv('d_family.csv', header=False, index=False)
    j_family_reindex.to_csv('j_family.csv', header=False, index=False)

    # # store dictionaries in a list and output to pickle file
    # gene_usage = [v_family, d_family, j_family]     ## excluding v_genes, d_genes, j_genes
    # with open('gene_usage_' + str(metadata[1] + '_' + str(metadata[2] + '_' + str(metadata[3]))) + '.pkl', 'wb') as f:
    #     pickle.dump(gene_usage, f)

calc_sample_stats(sample_meta, counts)