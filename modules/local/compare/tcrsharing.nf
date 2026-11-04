process TCRSHARING_CALC {
    label 'process_low'

    input:
    path concat_cdr3

    output:
    path "cdr3_sharing.tsv", emit: "shared_cdr3"
    path "sample_mapping.tsv", emit: "sample_mapping"

    script:
    """
    python - <<EOF
    import numpy as np
    import pandas as pd

    # Load data
    df = pd.read_csv("${concat_cdr3}", sep="\t")

    # Remove rows where pgen = 0
    df = df[df['pgen'] != 0]

    # Map sample to integer codes
    df['sample'] = df['sample'].astype('category')
    df['sample_id'] = df['sample'].cat.codes + 1

    # Export mapping (uses category lookup directly)
    sample_map_df = pd.DataFrame({
        'patient': df['sample'].cat.categories,
        'sample_id': np.arange(1, len(df['sample'].cat.categories) + 1)
    })
    sample_map_df.to_csv("sample_mapping.tsv", sep="\t", index=False)

    # Get unique sample_ids per CDR3b — vectorized
    grouped = (
        df.groupby('junction_aa')
        .agg(
            sample_id=('sample_id', 'unique'),
            pgen=('pgen', 'first'),
            log10_pgen=('log10_pgen', 'first')
        )
        .reset_index()
    )

    # Calculate counts
    grouped['total_samples'] = grouped['sample_id'].apply(len)
    grouped['samples_present'] = grouped['sample_id'].apply(
        lambda arr: ",".join(arr.astype(str))
    )

    # Drop raw list
    final_df = grouped[['junction_aa', 'pgen', 'log10_pgen', 'total_samples', 'samples_present']]
    final_df = final_df.sort_values(by="total_samples", ascending=False)

    # Export final list
    final_df.to_csv("cdr3_sharing.tsv", sep="\t", index=False)
    EOF
    """
}