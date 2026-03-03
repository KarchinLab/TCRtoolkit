#!/usr/bin/env python3

"""
gliph2_preprocess.py
Input: adaptive TSV files
Output: $concatenated_cdr3.txt
"""

# Import modules
import argparse
import pandas as pd

def main():
    # Initialize the parser
    parser = argparse.ArgumentParser(description="Take positional args")

    # Add positional arguments
    parser.add_argument("samplesheet")

    # Parse the arguments
    args = parser.parse_args()

    # Print the arguments
    print("samplesheet: ", args.samplesheet)
    samplesheet = pd.read_csv(args.samplesheet, header=0)
    dfs = []
    for _, row in samplesheet.iterrows():
        df = pd.read_csv(
            row['file'],
            sep="\t",
            usecols=[
                'amino_acid', #'junction_aa',
                'v_gene',
                'j_gene',
                'productive_frequency', #'duplicate_count',
                # 'productive'
            ]
        )
        df = df.rename(columns={
            'amino_acid': 'CDR3b',
            'v_gene': 'v_call',
            'j_gene': 'j_call',
            'productive_frequency': 'duplicate_frequency_percent'
        })

        # Retain only productive CDR3 sequences
        df = df[
            # (df['productive']) &
            (df['CDR3b'].notna()) # &
            # (df['v_call'].notna()) # also remove rows with a CDR3 sequence but no Vgene called
        ]

        df['sample'] = row['sample']
        df = df[['CDR3b', 'v_call', 'j_call', 'duplicate_frequency_percent', 'sample']]

        dfs.append(df)

    df_combined = pd.concat(dfs, ignore_index=True)

    # Rename columns as required
    # df_combined = df_combined.rename(columns={
    #     'junction_aa': 'CDR3b',
    #     'v_call': 'TRBV',
    #     'j_call': 'TRBJ',
    #     'duplicate_count': 'counts'
    # })

    df_combined.to_csv(f"concatenated_cdr3.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()