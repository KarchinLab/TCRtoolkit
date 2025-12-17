#!/usr/bin/env python3

"""
gliph2_preprocess.py
Input: adaptive TSV files
Output: $concatenated_cdr3.txt
"""

# Import modules
import argparse
import os
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
        # Read the TSV file into a dataframe
        file_path = str(row['file'])
        df = pd.read_csv(file_path, sep="\t", header=0)

        # Add patient column
        df['sample'] = row['sample']

        # Select relevant columns
        df = df[['junction_aa', 'v_call', 'j_call', 'duplicate_count', 'sample']]
        dfs.append(df)


    # Concatenate all the dataframes into one
    df_combined = pd.concat(dfs)

    # Rename columns as required
    df_combined = df_combined.rename(columns={
        'junction_aa': 'CDR3b',
        'v_call': 'TRBV',
        'j_call': 'TRBJ',
        'duplicate_count': 'counts'
    })
    df_combined = df_combined[df_combined['CDR3b'].notna()]

    df_combined.to_csv(f"concatenated_cdr3.txt", sep="\t", index=False, header=True)

if __name__ == "__main__":
    main()