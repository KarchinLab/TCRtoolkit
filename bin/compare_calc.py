#!/usr/bin/env python3
"""
Description: Calculate overlap measures between TCR repertoires
Author: Dylan Tamayo, Domenick Braccia
"""

import argparse
import pandas as pd
import numpy as np

# -------------------------
# Similarity functions
# -------------------------
def jaccard_index(set1, set2):
    union = len(set1 | set2)
    return len(set1 & set2) / union if union else 0.0


def sorensen_index(set1, set2):
    denom = len(set1) + len(set2)
    return (2 * len(set1 & set2) / denom) if denom else 0.0


def morisita_horn_index(counts1, counts2):
    X = counts1.sum()
    Y = counts2.sum()

    if X == 0 or Y == 0:
        return 0.0

    prod_sum = np.sum(counts1 * counts2)
    lambda1 = np.sum(counts1 ** 2) / (X ** 2)
    lambda2 = np.sum(counts2 ** 2) / (Y ** 2)

    return (2 * prod_sum) / ((lambda1 + lambda2) * X * Y)

if __name__ == "__main__":
    # -------------------------
    # Argument parsing
    # -------------------------
    parser = argparse.ArgumentParser(
        description="Calculate overlap metrics for TCR repertoires"
    )
    parser.add_argument(
        "-s", "--sample_utf8",
        required=True,
        help="Samplesheet CSV passed from Nextflow"
    )
    args = parser.parse_args()


    # -------------------------
    # Load samplesheet
    # -------------------------
    sample_df = pd.read_csv(args.sample_utf8)

    samples = sample_df["sample"].tolist()
    files = sample_df["file"].tolist()
    n = len(samples)

    print(f"Loaded {n} samples")

    # -------------------------
    # Preload data structures
    # -------------------------
    junction_sets = {}
    count_vectors = {}

    for sample, file in zip(samples, files):
        df = pd.read_csv(file, sep="\t", usecols=["junction_aa", "duplicate_count"])
        df = df.dropna(subset=["junction_aa"])

        # Set for presence/absence metrics
        junction_sets[sample] = set(df["junction_aa"])

        # Counts for Morisita–Horn
        count_vectors[sample] = (
            df.groupby("junction_aa")["duplicate_count"]
            .sum()
        )


    # -------------------------
    # Align count vectors across union space
    # -------------------------
    all_junctions = sorted(
        set().union(*junction_sets.values())
    )

    for sample in samples:
        count_vectors[sample] = (
            count_vectors[sample]
            .reindex(all_junctions, fill_value=0)
            .to_numpy()
        )


    # -------------------------
    # Initialize matrices
    # -------------------------
    jaccard_mat = np.zeros((n, n))
    sorensen_mat = np.zeros((n, n))
    morisita_mat = np.zeros((n, n))


    # -------------------------
    # Compute upper triangle only
    # -------------------------
    print("Calculating overlap metrics...")

    for i in range(n):
        s1 = samples[i]
        set1 = junction_sets[s1]
        counts1 = count_vectors[s1]

        # Diagonal
        jaccard_mat[i, i] = 1.0
        sorensen_mat[i, i] = 1.0
        morisita_mat[i, i] = 1.0

        for j in range(i + 1, n):
            s2 = samples[j]

            j_val = jaccard_index(set1, junction_sets[s2])
            s_val = sorensen_index(set1, junction_sets[s2])
            m_val = morisita_horn_index(counts1, count_vectors[s2])

            jaccard_mat[i, j] = jaccard_mat[j, i] = j_val
            sorensen_mat[i, j] = sorensen_mat[j, i] = s_val
            morisita_mat[i, j] = morisita_mat[j, i] = m_val


    # -------------------------
    # Write outputs
    # -------------------------
    index_names = samples

    pd.DataFrame(
        jaccard_mat, index=index_names, columns=index_names
    ).to_csv("jaccard_mat.csv")

    pd.DataFrame(
        sorensen_mat, index=index_names, columns=index_names
    ).to_csv("sorensen_mat.csv")

    pd.DataFrame(
        morisita_mat, index=index_names, columns=index_names
    ).to_csv("morisita_mat.csv")

    print("Finished writing all matrices")