#!/usr/bin/env python3
"""
Optimized TCR repertoire overlap (single-file input, vectorized)
"""

import argparse
import pandas as pd
import numpy as np


# -------------------------
# Morisita–Horn (vectorized)
# -------------------------
def morisita_horn_matrix(M):
    """
    M: (n_samples x n_features) count matrix
    Returns: (n x n) similarity matrix
    """
    X = M.sum(axis=1, keepdims=True)  # (n,1)

    # Avoid division by zero
    X_safe = np.where(X == 0, 1, X)

    lambda_vals = (M ** 2).sum(axis=1, keepdims=True) / (X_safe ** 2)

    # Pairwise dot products
    prod = M @ M.T  # (n x n)

    denom = (lambda_vals + lambda_vals.T) * (X @ X.T)

    # Avoid division by zero
    denom = np.where(denom == 0, 1, denom)

    return (2 * prod) / denom


# -------------------------
# Main
# -------------------------
if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Optimized overlap metrics for TCR repertoires"
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Concatenated TSV with sample, junction_aa, duplicate_count"
    )
    parser.add_argument(
        "-p", "--patient",
        required=True,
        help="String indicating patient of samples"
    )
    args = parser.parse_args()

    # -------------------------
    # Load + clean
    # -------------------------
    df = pd.read_csv(
        args.input,
        sep="\t",
        usecols=["sample", "junction_aa", "duplicate_count"]
    )

    df = df.dropna(subset=["sample", "junction_aa"])

    df["duplicate_count"] = pd.to_numeric(
        df["duplicate_count"], errors="coerce"
    ).fillna(0)

    # -------------------------
    # Pre-aggregate (critical optimization)
    # -------------------------
    df = (
        df.groupby(["sample", "junction_aa"], as_index=False)["duplicate_count"]
        .sum()
    )

    samples = df["sample"].unique()
    n = len(samples)

    print(f"Loaded {n} samples")

    # -------------------------
    # Build count matrix
    # -------------------------
    pivot = df.pivot(
        index="sample",
        columns="junction_aa",
        values="duplicate_count"
    ).fillna(0)

    # Ensure consistent ordering
    pivot = pivot.loc[samples]

    M = pivot.to_numpy(dtype=float)

    # -------------------------
    # Presence/absence matrix
    # -------------------------
    B = (M > 0).astype(np.int32)

    # -------------------------
    # Jaccard (vectorized)
    # -------------------------
    intersection = B @ B.T  # shared junction counts

    row_sums = B.sum(axis=1, keepdims=True)
    union = row_sums + row_sums.T - intersection

    union = np.where(union == 0, 1, union)

    jaccard = intersection / union

    # -------------------------
    # Sørensen (Dice)
    # -------------------------
    denom = row_sums + row_sums.T
    denom = np.where(denom == 0, 1, denom)

    sorensen = (2 * intersection) / denom

    # -------------------------
    # Morisita–Horn (vectorized)
    # -------------------------
    morisita = morisita_horn_matrix(M)

    # -------------------------
    # Fix diagonals
    # -------------------------
    np.fill_diagonal(jaccard, 1.0)
    np.fill_diagonal(sorensen, 1.0)
    np.fill_diagonal(morisita, 1.0)

    # -------------------------
    # Write outputs
    # -------------------------
    patient = args.patient
    pd.DataFrame(jaccard, index=samples, columns=samples).to_csv(f"{patient}/jaccard_mat.csv")
    pd.DataFrame(sorensen, index=samples, columns=samples).to_csv(f"{patient}/sorensen_mat.csv")
    pd.DataFrame(morisita, index=samples, columns=samples).to_csv(f"{patient}/morisita_mat.csv")

    print("Finished writing all matrices")