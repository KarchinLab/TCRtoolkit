#!/usr/bin/env python3

import argparse
import json
import logging

import pandas as pd

def configure_logging(args):
    """
    Configure logging based on the given arguments.
    """
    logger = logging.getLogger("pseudobulk")  # or use __name__ for module-level granularity

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


# Helper function for fields expected to have a single unique value (ignoring NaNs)
def assert_single_unique(x, field):
    non_na_values = set(x.dropna())
    if len(non_na_values) > 1:
        raise ValueError(f"Expected exactly one unique value for '{field}', but got: {non_na_values}")
    return next(iter(non_na_values)) if non_na_values else None

def airr_compliant_gene(gene):
    if pd.isna(gene):
        return None
    return gene + '*00' if '*' not in gene else gene

def normalize_bool(val):
    true_vals = {"true", "t", "1", "yes"}
    false_vals = {"false", "f", "0", "no"}
    
    if pd.isna(val):
        return pd.NA
    val_str = str(val).strip().lower()
    if val_str in true_vals:
        return "true"
    elif val_str in false_vals:
        return "false"
    return pd.NA 

def pseudobulk(input_df, basename, airr_schema):
    # Load data, filter out TRA calls
    df = input_df
    df_trb = df[df["v_call"].str.startswith("TRB") & df["j_call"].str.startswith("TRB")]

    # Define required aggregations
    agg_columns = {
        "cell_id": pd.Series.nunique,
        "junction": lambda x: assert_single_unique(x, "junction"),
        "junction_aa": lambda x: assert_single_unique(x, "junction_aa"),
        "v_call": lambda x: assert_single_unique(x, "v_call"),
        "j_call": lambda x: assert_single_unique(x, "j_call"),
        "consensus_count": "sum",
        "duplicate_count": "sum"
    }

    # Identify columns not already aggregated
    excluded_cols = set(agg_columns) | {"sequence"}
    remaining_cols = [col for col in df_trb.columns if col not in excluded_cols]

    # Group by sequence
    grouped = df_trb.groupby("sequence")

    # Determine which of the remaining columns have consistent values
    consensus_cols = []
    for col in remaining_cols:
        if (grouped[col].nunique(dropna=True) <= 1).all():
            consensus_cols.append(col)

    # Add those columns to the aggregation with the unique value
    for col in consensus_cols:
        agg_columns[col] = lambda x, col=col: assert_single_unique(x, col)

    # Aggregate
    bulk_df = grouped.agg(agg_columns).reset_index()

    # Rename and reorder columns
    bulk_df.rename(columns={"cell_id": "cell_count"}, inplace=True)
    bulk_df["sequence_id"] = [f"{basename}|sequence{i+1}" for i in range(len(bulk_df))]
    total_duplicate = bulk_df['duplicate_count'].sum()
    bulk_df['duplicate_frequency_percent'] = bulk_df['duplicate_count'] / total_duplicate * 100

    # Load schema
    schema = airr_schema

    required = schema["Rearrangement"]["required"]
    properties = list(schema["Rearrangement"]["properties"].keys())

    # Add missing required fields as blank
    for col in required:
        if col not in bulk_df.columns:
            bulk_df[col] = pd.NA

    # Apply to relevant columns if they exist
    for col in ['v_call', 'd_call', 'j_call', 'c_call']:
        if col in df.columns:
            bulk_df[col] = bulk_df[col].apply(airr_compliant_gene)

    # Identify all boolean columns from schema
    bool_cols = [
        name for name, spec in schema["Rearrangement"]["properties"].items()
        if spec.get("type") == "boolean"
    ]

    bool_cols = bool_cols + ['is_cell']
    for col in bool_cols:
        if col in bulk_df.columns:
            bulk_df[col] = bulk_df[col].apply(normalize_bool)

    # Order columns
    ordered_cols = [col for col in properties if col in bulk_df.columns]
    extra_cols = [col for col in bulk_df.columns if col not in ordered_cols]
    bulk_df = bulk_df[ordered_cols + extra_cols]

    output_path = f"{basename}_pseudobulk.tsv"
    bulk_df.to_csv(output_path, sep='\t', index=False)
    logger.info(f"Pseudobulked AIRR file written to: {output_path}")


def add_phenotype_info(input_df, basename, sobject_gex):
    import numpy as np
    import pandas as pd
    from rpy2.robjects import r
    from rpy2.robjects import pandas2ri
    import rpy2.robjects as ro
    from rpy2.robjects.conversion import localconverter

    # Read Seurat Object
    readRDS = r['readRDS']
    seurat_obj = readRDS(sobject_gex)
    # Extract and Convert Metadata into a pandas DataFrame
    with localconverter(ro.default_converter + pandas2ri.converter):
        cell_metadata_pd = ro.conversion.rpy2py(seurat_obj.slots['meta.data'])

    cell_metadata_pd["cell_id"] = cell_metadata_pd["library"] + "_" + cell_metadata_pd["cell_id"]
    cell_metadata_pd = cell_metadata_pd.loc[:,["cell_id","phenotype"]]
    cell_metadata_pd['phenotype'] = cell_metadata_pd['phenotype'].str.replace(' ', '_')

    # ------ Fix TCR Barcodes ------ #
    input_df["cell_id"] = basename + "_" + input_df["cell_id"]

    # ------ Merge ------ #
    merged_df = pd.merge(input_df, cell_metadata_pd, on="cell_id", how="left")
    return merged_df

def pseudobulk_phenotype(input_df, basename, airr_schema):
    # Load data, filter out TRA calls
    df = input_df
    df_trb = df[df["v_call"].str.startswith("TRB") & df["j_call"].str.startswith("TRB")]

    # Define required aggregations
    agg_columns = {
        "cell_id": pd.Series.nunique,
        "junction": lambda x: assert_single_unique(x, "junction"),
        "junction_aa": lambda x: assert_single_unique(x, "junction_aa"),
        "v_call": lambda x: assert_single_unique(x, "v_call"),
        "j_call": lambda x: assert_single_unique(x, "j_call"),
        "consensus_count": "sum",
        "duplicate_count": "sum"
    }

    # Identify columns not already aggregated
    excluded_cols = set(agg_columns) | {"sequence", "phenotype"}
    remaining_cols = [col for col in df_trb.columns if col not in excluded_cols]

    # Group by sequence
    grouped = df_trb.groupby(["sequence", "phenotype"])

    # Determine which of the remaining columns have consistent values
    consensus_cols = []
    for col in remaining_cols:
        if (grouped[col].nunique(dropna=True) <= 1).all():
            consensus_cols.append(col)

    # Add those columns to the aggregation with the unique value
    for col in consensus_cols:
        agg_columns[col] = lambda x, col=col: assert_single_unique(x, col)

    # Aggregate
    bulk_df = grouped.agg(agg_columns).reset_index()

    # Rename and reorder columns
    bulk_df.rename(columns={"cell_id": "cell_count"}, inplace=True)
    bulk_df["sequence_id"] = [f"{basename}|sequence{i+1}|{bulk_df['phenotype'].iloc[i]}" for i in range(len(bulk_df))]
    total_duplicate = bulk_df['duplicate_count'].sum()
    bulk_df['duplicate_frequency_percent'] = bulk_df['duplicate_count'] / total_duplicate * 100

    # Load schema
    schema = airr_schema

    required = schema["Rearrangement"]["required"]
    properties = list(schema["Rearrangement"]["properties"].keys())

    # Add missing required fields as blank
    for col in required:
        if col not in bulk_df.columns:
            bulk_df[col] = pd.NA

    # Apply to relevant columns if they exist
    for col in ['v_call', 'd_call', 'j_call', 'c_call']:
        if col in df.columns:
            bulk_df[col] = bulk_df[col].apply(airr_compliant_gene)

    # Identify all boolean columns from schema
    bool_cols = [
        name for name, spec in schema["Rearrangement"]["properties"].items()
        if spec.get("type") == "boolean"
    ]

    bool_cols = bool_cols + ['is_cell']
    for col in bool_cols:
        if col in bulk_df.columns:
            bulk_df[col] = bulk_df[col].apply(normalize_bool)

    # Order columns
    ordered_cols = [col for col in properties if col in bulk_df.columns]
    extra_cols = [col for col in bulk_df.columns if col not in ordered_cols]
    bulk_df = bulk_df[ordered_cols + extra_cols]

    # Iterate through each unique phenotype and save to a separate file
    for phenotype in bulk_df['phenotype'].unique():
        # Filter the dataframe for the current phenotype
        phenotype_df = bulk_df[bulk_df['phenotype'] == phenotype]

        # Define the output path based on basename and phenotype
        output_path = f"{basename}_{phenotype}_pseudobulk_phenotype.tsv"

        # Save the filtered dataframe to a TSV file
        phenotype_df.to_csv(output_path, sep='\t', index=False)
        logger.info(f"Pseudobulked AIRR file for phenotype '{phenotype}' written to: {output_path}")


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
        "airr_schema",
        type=str,
        help="Path to airr json schema for validating adaptive conversion.",
    )

    parser.add_argument(
        "--phenotype",
        action="store_true",
        help="A boolean flag to indicate wether the pseudobulk approach will consider phenotype.",
    )

    parser.add_argument(
        "--sobject_gex",
        type=str,
        help="Path to Seurat Object, 'phenotype' should be a column in meta.data slot.",
    )

    add_logging_args(parser)
    args = parser.parse_args()

    configure_logging(args)
    logger = logging.getLogger("convert")

    logger.info("Beginning Adaptive Immunoseq to AIRR format conversion...")
    if args.verbose:
        logger.debug("Verbose logging enabled")

    logger.info(f"sample_tsv: {args.sample_tsv}")
    logger.info(f"basename: {args.basename}")
    logger.info(f"airr_schema: {args.airr_schema}")
    
    input_df = pd.read_csv(args.sample_tsv, sep="\t")
    with open(args.airr_schema) as f:
        airr_schema = json.load(f)
    basename = args.basename

    pseudobulk(input_df, basename, airr_schema)

    # Pseudobulk by phenotype (if available)
    if args.phenotype:
        input_df = add_phenotype_info(input_df, basename, args.sobject_gex)  # Add single-cell phenotype info from Seurat Object.
        pseudobulk_phenotype(input_df, basename, airr_schema)
    