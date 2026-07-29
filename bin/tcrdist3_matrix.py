#!/usr/bin/env python3

import argparse
import logging
import re

import h5py
import numpy as np
import pandas as pd
import Levenshtein
import scipy.sparse as sp
from tcrdist.repertoire import TCRrep

def configure_logging(args):
    """
    Configure logging based on the given arguments.
    """
    logger = logging.getLogger("tcrdist")  # or use __name__ for module-level granularity

    ch = logging.StreamHandler()

    if args.verbose:
        logger.setLevel(logging.DEBUG)
        ch.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)
        ch.setLevel(logging.INFO)

    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    ch.setFormatter(formatter)
    if not logger.hasHandlers():
        logger.addHandler(ch)

def add_logging_args(parser):
    parser.add_argument(
        "-v",
        "--verbose",
        default=False,
        action="store_true",
        help="Enable verbose logging",
    )

def transform_trbv(trbv):
    """Convert gene names from Adaptive ImmunoSEQ to IMGT format."""
    if not isinstance(trbv, str):
        return trbv  # Return as-is if not a string

    # Convert locus name
    trbv = trbv.replace("TCRBV", "TRBV")

    # Remove zero padding from gene name (TCRBV07 to TRBV7)
    trbv = re.sub(r'(?<=TRBV)0*(\d+)', r'\1', trbv)  

    # Remove zero padding from subgroup (TCRBV7-02 to TRBV7-2)
    trbv = re.sub(r'-(0\d+)', lambda m: f'-{int(m.group(1))}', trbv)  

    # Convert "-orXX_XX" to "/OR#-#" for orphon genes
    trbv = re.sub(r'-or0?(\d+)_0?(\d+)', r'/OR\1-\2', trbv)

    # Add *01 if allele group not specified
    if not re.search(r'\*\d{2}$', trbv):
        trbv += "*01"

    return trbv

def remove_locus(gene_name):
    """Remove the -## gene position from TRBV names unless the gene contains /OR."""
    if '/OR' in gene_name:
        return gene_name
    else:
        return re.sub(r'-(\d+)\*', '*', gene_name)

def split_and_check_genes(gene_name):
    """Split combined TCRBV genes (e.g., TCRBV06-02/06-03*01 to TCRBV06-02*01 and TCRBV06-03*01."""
    if '/' in gene_name and not re.search(r'/OR\d+-\d+', gene_name):  # Ensure it's not an orphon
        base, allele = gene_name.split("*") if "*" in gene_name else (gene_name, "01")
        prefix_match = re.match(r"(TCRBV\d+)", gene_name)
        prefix = prefix_match.group(1) if prefix_match else "TCRBV"  # Fallback just in case
        genes = base.split("/")  # Split the genes
        return [f"{prefix}-{g.split('-')[-1]}*{allele}" for g in genes]
    return [gene_name]


def find_matching_gene(row, db):
    # Collect all possible genes from vMaxResolved and vGeneNameTies
    possible_genes = set()  # Use a set to avoid duplicates

    if pd.notna(row["vMaxResolved"]):
        possible_genes.add(row["vMaxResolved"])  # Always include vMaxResolved

    if pd.notna(row["vGeneNameTies"]):
        possible_genes.update(row["vGeneNameTies"].split(","))  # Add vGeneNameTies genes

    for gene in possible_genes:
        # If the gene contains multiple variants (e.g., TCRBV03-01/03-02*01), split and check both
        if "/" in gene and not re.search(r"/OR\d+-\d+", gene):  # Avoid /OR cases
            sub_genes = split_and_check_genes(gene)
            for sub_gene in sub_genes:
                sub_gene = transform_trbv(sub_gene)  # Ensure correct *0# format
                if sub_gene in db["id"].values:
                    return sub_gene

        # Direct match in db
        transform_gene = transform_trbv(gene)
        if transform_gene in db["id"].values:
            return transform_gene

        # Try removing -## and checking again
        modified_gene = remove_locus(transform_gene)
        if modified_gene in db["id"].values:
            return modified_gene

    transform_row = transform_trbv(row["vMaxResolved"])
    print(f'No match found for {transform_row}')

    return transform_row  # Return original vMaxResolved if no match is found


if __name__ == "__main__":
    # Parse input arguments
    parser = argparse.ArgumentParser(description="Take positional args")

    parser.add_argument(
        "sample_tsv",
        type=str,
        help="Path to TSV file containing TCR count table.",
    )

    parser.add_argument(
        "basename",
        type=str,
        help="Sample ID name for output files.",
    )

    parser.add_argument(
        "matrix_sparsity",
        type=str,
        choices=["sparse", "full"],
        help="Matrix type to compute: 'sparse' for memory-efficient npz format or 'full' for complete CSV matrix.",
    )

    parser.add_argument(
        "distance_metric",
        type=str,
        choices=["levenshtein", "tcrdist"],
        help="Distance metric to use for computing pairwise distances between TCRs.",
    )

    parser.add_argument(
        "ref_database",
        type=str,
        help="Path to reference database file for matching genes from input file.",
    )

    parser.add_argument(
        "cores",
        type=int,
        help="Number of CPU cores to use for parallel computation.",
    )
    add_logging_args(parser)
    args = parser.parse_args()

    configure_logging(args)
    logger = logging.getLogger("tcrdist")

    logger.info("Beginning tcrdist3 matrix generation...")
    if args.verbose:
        logger.debug("Verbose logging enabled")

    logger.info(f"sample_tsv: {args.sample_tsv}")
    logger.info(f"basename: {args.basename}")
    logger.info(f"matrix_sparsity: {args.matrix_sparsity}")
    logger.info(f"distance_metric: {args.distance_metric}")
    logger.info(f"ref_database: {args.ref_database}")
    logger.info(f"cores: {args.cores}")

    basename = args.basename

    # --- 1. Convert Adaptive output to tcrdist db format ---
    db = pd.read_table(args.ref_database, delimiter = '\t')

    db = db[db['organism']=='human']

    df = pd.read_table(args.sample_tsv, delimiter = '\t')

    df = df[['sequence', 'junction_aa', 'v_call', 'duplicate_count']]

    df = df.rename(columns={'sequence': 'cdr3_b_nucseq',
                        'junction_aa': 'cdr3_b_aa',
                        'v_call': 'v_b_gene',
                        'duplicate_count': 'count'
                        })

    df = df[df['cdr3_b_aa'].notna()]
    df = df[df['v_b_gene'].notna()]

    # If allele information is not specified (as indicated by *00), replace with *01
    df['v_b_gene'] = df['v_b_gene'].apply(lambda x: x.replace('*00', '*01'))

    # Single-cell pseudobulk carries no CDR3 nucleotide sequence. tcrdist3 uses cdr3_b_nucseq
    # as a clone grouping/index column, so an all-empty nucseq would drop every clone
    # ("N of N were not captured" -> empty matrix). Beta-chain tcrdist distances are computed
    # from cdr3_b_aa and V-gene-derived CDR1/CDR2, not the nucleotide sequence, so drop the
    # column when it is empty. (Bulk / AIRR input has real nucseq and is unaffected.)
    if 'cdr3_b_nucseq' in df.columns:
        _ns = df['cdr3_b_nucseq'].astype('string').fillna('').str.strip()
        if (_ns == '').all():
            df = df.drop(columns=['cdr3_b_nucseq'])

    # --- 2. Calculate distance matrix ---
    # Levenshtein distance matrix
    if args.distance_metric == "levenshtein":
        def my_own_metric(s1,s2):   
            return Levenshtein.distance(s1,s2)

        tr = TCRrep(cell_df = df,
                    organism = 'human',
                    chains = ['beta'],
                    db_file = 'alphabeta_gammadelta_db.tsv',
                    use_defaults=False,
                    compute_distances = False)

        metrics_b = {
            "cdr3_b_aa" : my_own_metric,
            "pmhc_b_aa" : my_own_metric,
            "cdr2_b_aa" : my_own_metric,
            "cdr1_b_aa" : my_own_metric }

        weights_b = { 
            "cdr3_b_aa" : 1,
            "pmhc_b_aa" : 1,
            "cdr2_b_aa" : 1,
            "cdr1_b_aa" : 1}

        kargs_b = {  
            'cdr3_b_aa' : 
                {'use_numba': False},
            'pmhc_b_aa' : {
                'use_numba': False},
            'cdr2_b_aa' : {
                'use_numba': False},
            'cdr1_b_aa' : {
                'use_numba': False}
            }

        tr.metrics_b = metrics_b
        tr.weights_b = weights_b
        tr.kargs_b = kargs_b
        radius = 6

    # Default tcrdist3 distance matrix
    else:
        tr = TCRrep(cell_df = df,
                    organism = 'human',
                    chains = ['beta'],
                    db_file = 'alphabeta_gammadelta_db.tsv',
                    compute_distances = False)
        radius = 50

    tr.cpus = args.cores
    clone_df = tr.clone_df
    clone_df.to_csv(f"{basename}_clone_df.csv", index=False)

    # Full matrix
    if args.matrix_sparsity == "full":
        tr.compute_distances()
        max_val = tr.pw_beta.max()
        np.savetxt(f"{basename}_distance_matrix.csv", tr.pw_beta, delimiter=",", fmt="%d")
    # Sparse matrix
    else:
        tr.compute_sparse_rect_distances(radius = radius, chunk_size = 100)
        max_val = tr.rw_beta.max()
        with h5py.File(f"{basename}_distance_matrix.hdf5", "w") as f:
            f.create_dataset("data", data=tr.rw_beta.data)
            f.create_dataset("indices", data=tr.rw_beta.indices)
            f.create_dataset("indptr", data=tr.rw_beta.indptr)
            f.create_dataset("shape", data=tr.rw_beta.shape)

    with open(f"matrix_maximum_value.txt", "w") as f:
        f.write(str(max_val))