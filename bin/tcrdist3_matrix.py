#!/usr/bin/env python3

import argparse
import logging
import re
import resource
import sys
import time

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


def peak_rss_mb():
    """Peak resident memory of this process so far, in MB (Linux reports KB)."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def cluster_map(clone_df, mat, radius, patient):
    """
    TCRdist clusters at `radius`, written once here so downstream steps (CLUSTER_TO_SC,
    CLUSTER_ROLLUP) read a small table instead of every patient's raw distance matrix.

    Same definition as enrich_seurat.R had: clone rows within `radius` are linked, rows
    sharing a CDR3b are tied together (cells are matched on CDR3b), and every connected
    component with >= 2 distinct CDR3b is a cluster. In the sparse layout an absent entry
    is beyond the cutoff and a true zero distance is stored as -1. The matrix stays sparse.
    """
    from scipy.sparse.csgraph import connected_components

    n = len(clone_df)
    if sp.issparse(mat):
        coo = mat.tocoo()
        d = coo.data.copy()
        d[d == -1] = 0
        keep = (d <= radius) & (coo.row < coo.col)
        ei, ej = coo.row[keep], coo.col[keep]
    else:
        ei, ej = np.nonzero(np.triu(mat <= radius, k=1))
    cdr3 = clone_df["cdr3_b_aa"].astype(str).str.upper().to_numpy()
    cdr3_codes, cdr3_uniq = pd.factorize(cdr3)
    # Nodes 0..n-1 are clone rows, n.. are CDR3b tie nodes.
    rows = np.concatenate([ei, np.arange(n)])
    cols = np.concatenate([ej, n + cdr3_codes])
    adj = sp.coo_matrix((np.ones(len(rows), dtype=np.int8), (rows, cols)),
                        shape=(n + len(cdr3_uniq),) * 2)
    _, labels = connected_components(adj, directed=False)
    comp = labels[:n]

    out = pd.DataFrame({
        "patient": patient,
        "clone_index": np.arange(n),
        "sample": clone_df["sample"].to_numpy() if "sample" in clone_df else "",
        "cdr3_b_aa": cdr3,
        "v_b_gene": clone_df["v_b_gene"].to_numpy() if "v_b_gene" in clone_df else "",
        "count": clone_df["count"].to_numpy() if "count" in clone_df else np.nan,
        "_comp": comp,
    })
    n_cdr3 = out.groupby("_comp")["cdr3_b_aa"].nunique()
    n_rows = out.groupby("_comp").size()
    n_edges = pd.Series(comp[ei]).value_counts() if len(ei) else pd.Series(dtype=int)
    clustered = n_cdr3.index[n_cdr3 >= 2]
    # Labels numbered by first appearance, so they are stable for a given input order.
    order = pd.unique(comp[np.isin(comp, clustered)])
    label = {c: f"{patient}_T{i + 1}" for i, c in enumerate(order)}
    out["tcrdist_cluster"] = out["_comp"].map(label).fillna("")
    out["cluster_n_cdr3"] = out["_comp"].map(n_cdr3).where(out["tcrdist_cluster"] != "", 0).astype(int)
    out["cluster_n_clones"] = out["_comp"].map(n_rows).where(out["tcrdist_cluster"] != "", 0).astype(int)
    out["cluster_n_edges"] = out["_comp"].map(n_edges).fillna(0).where(out["tcrdist_cluster"] != "", 0).astype(int)
    return out.drop(columns="_comp"), int(len(ei))


def write_matrix(mat, prefix, sparsity):
    """Write a distance matrix in the layout the report notebooks and enrich_seurat.R read."""
    if sparsity == "full":
        np.savetxt(f"{prefix}.csv", mat, delimiter=",", fmt="%d")
    else:
        mat = sp.csr_matrix(mat)
        with h5py.File(f"{prefix}.hdf5", "w") as f:
            f.create_dataset("data", data=mat.data)
            f.create_dataset("indices", data=mat.indices)
            f.create_dataset("indptr", data=mat.indptr)
            f.create_dataset("shape", data=mat.shape)


if __name__ == "__main__":
    # Parse input arguments
    parser = argparse.ArgumentParser(
        description="TCRdist3 distances for one patient: every clone of every sample of "
                    "that patient against every other, so similar TCRs are found across "
                    "the patient's timepoints and tissues. Each sample's block of the "
                    "patient matrix is also written under the per-sample file names the "
                    "per-sample report sections read."
    )

    parser.add_argument(
        "basename",
        type=str,
        help="Patient ID used to name the patient-level output files.",
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

    parser.add_argument(
        "--radius",
        type=float,
        default=None,
        help="Distance cutoff: sparse matrices store only pairs within it, and the cluster map "
             "links clones within it (--tcrdist_radius). Default: 50 (tcrdist) / 6 (levenshtein), "
             "the previous hard-coded cutoffs.",
    )

    parser.add_argument(
        "--max_clones",
        type=int,
        default=0,
        help="Fail before any distance is computed when the patient has more clones than this "
             "(--tcrdist_max_patient_clones). 0 = no limit.",
    )

    parser.add_argument(
        "--input",
        nargs="+",
        required=True,
        metavar="SAMPLE=TSV",
        help="Per-sample TCR count tables of this patient, as sample=path pairs.",
    )
    add_logging_args(parser)
    args = parser.parse_args()

    configure_logging(args)
    logger = logging.getLogger("tcrdist")
    t_start = time.time()

    logger.info("Beginning tcrdist3 matrix generation...")
    if args.verbose:
        logger.debug("Verbose logging enabled")

    logger.info(f"patient: {args.basename}")
    logger.info(f"inputs: {args.input}")
    logger.info(f"matrix_sparsity: {args.matrix_sparsity}")
    logger.info(f"distance_metric: {args.distance_metric}")
    logger.info(f"ref_database: {args.ref_database}")
    logger.info(f"cores: {args.cores}")

    basename = args.basename

    # --- 1. Convert Adaptive output to tcrdist db format ---
    db = pd.read_table(args.ref_database, delimiter = '\t')

    db = db[db['organism']=='human']

    frames = []
    for item in args.input:
        sample, _, path = item.partition("=")
        part = pd.read_table(path, delimiter = '\t')
        part = part[['sequence', 'junction_aa', 'v_call', 'duplicate_count']]
        # 'sample' is an index column for TCRrep (every column except count / cell_id /
        # clone_id is), so the same TCR in two samples stays two clones.
        part['sample'] = sample
        frames.append(part)
    df = pd.concat(frames, ignore_index=True)
    samples = [item.partition("=")[0] for item in args.input]
    n_input_rows = int(len(df))

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
        radius = 6 if args.radius is None else args.radius

    # Default tcrdist3 distance matrix
    else:
        tr = TCRrep(cell_df = df,
                    organism = 'human',
                    chains = ['beta'],
                    db_file = 'alphabeta_gammadelta_db.tsv',
                    compute_distances = False)
        radius = 50 if args.radius is None else args.radius

    tr.cpus = args.cores
    clone_df = tr.clone_df.reset_index(drop=True)
    n_patient_clones = int(len(clone_df))
    n_unique_cdr3 = int(clone_df["cdr3_b_aa"].nunique())
    logger.info(f"n_input_rows={n_input_rows} n_unique_cdr3_beta={n_unique_cdr3} "
                f"n_patient_clones={n_patient_clones} radius={radius} peak_rss_mb={peak_rss_mb():.0f}")

    # Fail before building the matrix: memory grows with the number of clone pairs, so a
    # patient over the limit would otherwise die after hours of distance computation.
    if args.max_clones and n_patient_clones > args.max_clones:
        logger.error(f"TCRDIST_LIMIT patient={basename} n_patient_clones={n_patient_clones} "
                     f"exceeds --tcrdist_max_patient_clones={args.max_clones}; no matrix computed. "
                     f"Raise the limit (and memory) or downsample this patient.")
        sys.exit(3)

    clone_df.to_csv(f"{basename}_tcrdist_patient_clones.csv", index=False)

    # Full matrix
    if args.matrix_sparsity == "full":
        tr.compute_distances()
        mat = tr.pw_beta
        max_val = mat.max()
    # Sparse matrix (true zero distances are stored as -1, as tcrdist3 does)
    else:
        tr.compute_sparse_rect_distances(radius = radius, chunk_size = 100)
        mat = tr.rw_beta.tocsr()
        max_val = mat.max()
    write_matrix(mat, f"{basename}_tcrdist_patient_matrix", args.matrix_sparsity)
    nnz = int(mat.nnz) if sp.issparse(mat) else int(np.count_nonzero(mat <= radius))

    cmap, n_edges = cluster_map(clone_df, mat, radius, basename)
    cmap.to_csv(f"{basename}_tcrdist_cluster_map.tsv", sep="\t", index=False)
    logger.info(f"cluster map: {cmap['tcrdist_cluster'].replace('', np.nan).nunique()} clusters, "
                f"{n_edges} edges within radius {radius}")

    # --- 3. Per-sample blocks of the patient matrix ---
    # A distance depends only on the two TCRs compared, so a sample's block is exactly the
    # matrix a per-sample run would have produced.
    for sample in samples:
        idx = np.flatnonzero(clone_df['sample'].to_numpy() == sample)
        sample_clones = clone_df.iloc[idx].drop(columns=['sample'])
        sample_clones.to_csv(f"{sample}_clone_df.csv", index=False)
        block = mat[idx][:, idx]
        write_matrix(block, f"{sample}_distance_matrix", args.matrix_sparsity)

    with open(f"matrix_maximum_value.txt", "w") as f:
        f.write(str(max_val))

    # One row per patient; the workflow collects them into tcrdist_metrics.tsv.
    pd.DataFrame([{
        "patient": basename,
        "n_samples": len(samples),
        "n_input_rows": n_input_rows,
        "n_unique_cdr3_beta": n_unique_cdr3,
        "n_patient_clones": n_patient_clones,
        "matrix_sparsity": args.matrix_sparsity,
        "radius": radius,
        "nnz": nnz,
        "n_edges_within_radius": n_edges,
        "n_clusters": int(cmap["tcrdist_cluster"].replace("", np.nan).nunique()),
        "cpus": args.cores,
        "wall_time_s": round(time.time() - t_start, 1),
        "peak_rss_mb": round(peak_rss_mb(), 1),
    }]).to_csv(f"{basename}_tcrdist_metrics.tsv", sep="\t", index=False)
