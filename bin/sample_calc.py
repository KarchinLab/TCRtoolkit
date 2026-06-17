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
import re

def extract_trb_family(allele):
    if pd.isna(allele):
        return None
    match = re.match(r'(TRB[V|D|J])(\d+)', allele)
    return f"{match.group(1)}{match.group(2)}" if match else None

def calc_gene_family(sample_name, counts, gene_column, output_file):
    (
        counts[gene_column]
        .apply(extract_trb_family)
        .dropna()
        .value_counts()
        .rename_axis('gene')
        .reset_index(name='count')
        .assign(sample=sample_name)[['sample', 'gene', 'count']]
        .to_csv(output_file, index=False)
    )

def calc_sample_stats(sample_name, counts, pre_stats, output_file):
    """Calculate sample level statistics of TCR repertoire."""

    ## first pass stats
    clone_counts = counts['duplicate_count']
    clone_entropy = entropy(clone_counts, base=2)
    num_clones = len(clone_counts)
    num_TCRs = sum(clone_counts)
    clonality = 1 - clone_entropy / np.log2(num_clones)
    simpson_index = sum(clone_counts**2)/(num_TCRs**2)
    simpson_index_corrected = sum(clone_counts*(clone_counts-1))/(num_TCRs*(num_TCRs-1))

    # productive stats from pre-filter sidecar (input is productive-only)
    num_prod = int(pre_stats['productive_clones'])
    num_nonprod = int(pre_stats['nonproductive_clones'])
    total_pre_filter = int(pre_stats['total_clones'])
    pct_prod = num_prod / total_pre_filter if total_pre_filter > 0 else 0.0
    pct_nonprod = num_nonprod / total_pre_filter if total_pre_filter > 0 else 0.0

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

    # Add sample column
    df_stats.insert(0, 'sample', sample_name)

    # Save to CSV
    df_stats.to_csv(output_file, header=True, index=False)


def main():
    # initialize parser
    parser = argparse.ArgumentParser(description='Calculate clonality of a TCR repertoire')

    # add arguments
    parser.add_argument('-s', '--sample_name', 
                        metavar='sample_name', 
                        type=str, 
                        help='sample name')
    parser.add_argument('-c', '--count_table', 
                        metavar='count_table', 
                        type=str, 
                        help='counts file in TSV format')
    parser.add_argument('-p', '--pre_filter_stats',
                        metavar='pre_filter_stats',
                        type=str,
                        help='pre-filter stats CSV from ANNOTATE_PROCESS')

    args = parser.parse_args() 

    sample = args.sample_name
    counts = pd.read_csv(args.count_table, sep='\t')
    pre_filter_stats = pd.read_csv(args.pre_filter_stats).iloc[0]

    calc_gene_family(sample, counts, 'v_call', f'v_family_{sample}.csv')
    calc_gene_family(sample, counts, 'd_call', f'd_family_{sample}.csv')
    calc_gene_family(sample, counts, 'j_call', f'j_family_{sample}.csv')

    calc_sample_stats(sample, counts, pre_filter_stats, f'sample_stats_{sample}.csv')

if __name__ == "__main__":
    main()
