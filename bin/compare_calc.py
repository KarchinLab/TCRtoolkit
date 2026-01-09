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

    # Basic hygiene: drop rows missing sample or file
    sample_df = sample_df.dropna(subset=["sample", "file"])

    samples = sample_df["sample"].tolist()
    files = sample_df["file"].tolist()
    n = len(samples)

    print(f"Loaded {n} samples")

    # -------------------------
    # Preload data structures
    # -------------------------
    junction_sets = {}
    count_series = {}

    for sample, file in zip(samples, files):
        df = pd.read_csv(file, sep="\t", usecols=["junction_aa", "duplicate_count"])
        df = df.dropna(subset=["junction_aa"])

        # Ensure counts are numeric
        df["duplicate_count"] = pd.to_numeric(df["duplicate_count"], errors="coerce").fillna(0)

        # Set for presence/absence metrics
        junction_sets[sample] = set(df["junction_aa"])

        # Counts for Morisita–Horn as a pandas Series (index = junction_aa)
        counts = (
            df.groupby("junction_aa")["duplicate_count"]
            .sum()
        )
        # Ensure we have a Series with a named index
        if not isinstance(counts, pd.Series):
            counts = pd.Series(counts)
        count_series[sample] = counts

    # -------------------------
    # Build union index and align counts
    # -------------------------
    union_set = set().union(*junction_sets.values()) if junction_sets else set()
    all_junctions = pd.Index(sorted(union_set))

    aligned_vectors = {}
    for sample in samples:
        aligned = count_series[sample].reindex(all_junctions, fill_value=0)
        # Store as numpy for MH computation
        aligned_vectors[sample] = aligned.to_numpy(dtype=float)

    # -------------------------
    # Initialize matrices
    # -------------------------
    jaccard_mat = np.zeros((n, n), dtype=float)
    sorensen_mat = np.zeros((n, n), dtype=float)
    morisita_mat = np.zeros((n, n), dtype=float)

    # -------------------------
    # Compute upper triangle only
    # -------------------------
    print("Calculating overlap metrics...")

    for i in range(n):
        s1 = samples[i]
        set1 = junction_sets[s1]
        counts1 = aligned_vectors[s1]

        # Diagonal
        jaccard_mat[i, i] = 1.0
        sorensen_mat[i, i] = 1.0
        morisita_mat[i, i] = 1.0

        for j in range(i + 1, n):
            s2 = samples[j]
            j_val = jaccard_index(set1, junction_sets[s2])
            s_val = sorensen_index(set1, junction_sets[s2])
            m_val = morisita_horn_index(counts1, aligned_vectors[s2])

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