#!/usr/bin/env python3
"""
Description: this script calculates the clonality of a TCR repertoire

@author: Dylan Tamayo, Domenick Braccia
@contributor: elhanaty
"""

## import packages
import argparse
import pandas as pd
import numpy as np
from scipy.stats import entropy
import numpy as np
import csv
import re
import json

def extract_trb_family(allele):
    if pd.isna(allele):
        return None
    match = re.match(r'(TRB[V|D|J])(\d+)', allele)
    return f"{match.group(1)}{match.group(2)}" if match else None

def compute_gene_family_table(counts, col_name, all_families, sample_meta):
    fam_col = f"{col_name}FamilyName"
    counts[fam_col] = counts[col_name].apply(extract_trb_family)
    fam_df = counts[fam_col].value_counts(dropna=False).to_frame().T.sort_index(axis=1)
    fam_df = fam_df.reindex(columns=all_families, fill_value=0)

    for col in ['origin', 'timepoint', 'subject_id']:
        fam_df.insert(0, col, sample_meta[col])

    return fam_df

def calc_gene_family(counts, gene_column, family_prefix, max_index, output_file, meta_df):
    # Build list of all possible family names
    all_fams = [f'{family_prefix}{i}' for i in range(1, max_index + 1)]

    # Count usage
    fam_df = counts[gene_column].apply(extract_trb_family).value_counts(dropna=False).to_frame().T

    # Reindex to include all families
    fam_df = pd.DataFrame([fam_df.reindex(columns=all_fams, fill_value=0).iloc[0]]).reset_index(drop=True)

    # Add metadata columns
    fam_df = pd.concat([meta_df, fam_df], axis=1)

    fam_df.to_csv(output_file, index=False, header=False)

def calc_sample_stats(meta_df, counts):
    """Calculate sample level statistics of TCR repertoire."""

    ## first pass stats
    clone_counts = counts['duplicate_count']
    clone_entropy = entropy(clone_counts, base=2)
    num_clones = len(clone_counts)
    num_TCRs = sum(clone_counts)
    clonality = 1 - clone_entropy / np.log2(num_clones)
    simpson_index = sum(clone_counts**2)/(num_TCRs**2)
    simpson_index_corrected = sum(clone_counts*(clone_counts-1))/(num_TCRs*(num_TCRs-1))

    # count number of productive clones
    num_prod = sum(counts['productive'])
    num_nonprod = num_clones - num_prod
    pct_prod = num_prod / num_clones
    pct_nonprod = num_nonprod / num_clones

    ## cdr3 info
    cdr3_lens = counts['junction_aa_length']
    productive_cdr3_avg_len = np.mean([x*3 for x in cdr3_lens if x > 0])

    ## Calculate convergence for each T cell receptor
    aas = counts[counts.junction_aa.notnull()].junction_aa.unique()
    dict_df = {}
    for aa in aas:
        dict_df[aa] = {'counts': counts[counts.junction_aa == aa]}
        # append key value pair to dict_df[aa] with key convergence equal to the number of rows in counts
        dict_df[aa]['convergence'] = len(counts[counts.junction_aa == aa])

    ## calculate the number of covergent TCRs for each sample
    num_convergent = 0
    for aa in aas:
        if dict_df[aa]['convergence'] > 1:
            num_convergent += 1    

    ## calculate ratio of convergent TCRs to total TCRs
    ratio_convergent = num_convergent/len(aas)

    row_data = {
        'num_clones': num_clones,
        'num_TCRs': num_TCRs,
        'simpson_index': simpson_index,
        'simpson_index_corrected': simpson_index_corrected,
        'clonality': clonality,
        'num_prod': num_prod,
        'num_nonprod': num_nonprod,
        'pct_prod': pct_prod,
        'pct_nonprod': pct_nonprod,
        'productive_cdr3_avg_len': productive_cdr3_avg_len,
        'num_convergent': num_convergent,
        'ratio_convergent': ratio_convergent
    }

    # Convert to single-row dataframe
    df_stats = pd.DataFrame([row_data])

    # Add metadata columns
    df_stats = pd.concat([meta_df, df_stats], axis=1)

    # Save to CSV
    df_stats.to_csv('sample_stats.csv', header=False, index=False)


def main():
    # initialize parser
    parser = argparse.ArgumentParser(description='Calculate clonality of a TCR repertoire')

    # add arguments
    parser.add_argument('-s', '--sample_meta', 
                        metavar='sample_meta', 
                        type=str, 
                        help='sample metadata passed in as json format')
    parser.add_argument('-c', '--count_table', 
                        metavar='count_table', 
                        type=argparse.FileType('r'), 
                        help='counts file in TSV format')

    args = parser.parse_args() 

    ## convert metadata to list
    sample_meta = json.loads(args.sample_meta)

    # Read in the counts file
    counts = pd.read_csv(args.count_table, sep='\t', header=0)

    # Build metadata row from selected keys
    meta_keys = ['subject_id', 'timepoint', 'origin']
    meta_row = {k: sample_meta[k] for k in meta_keys}
    meta_df = pd.DataFrame([meta_row])

    calc_gene_family(counts, 'v_call', 'TRBV', 30, 'v_family.csv', meta_df)
    calc_gene_family(counts, 'd_call', 'TRBD', 2, 'd_family.csv', meta_df)
    calc_gene_family(counts, 'j_call', 'TRBJ', 2, 'j_family.csv', meta_df)

    # Build metadata row from selected keys
    meta_keys = ['sample', 'subject_id', 'timepoint', 'origin']
    meta_row = {k: sample_meta[k] for k in meta_keys}
    meta_df = pd.DataFrame([meta_row])
    
    calc_sample_stats(meta_df, counts)

if __name__ == "__main__":
    main()
