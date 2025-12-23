process TCRSHARING_CALC {
    label 'process_low'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    input:
    path concat_cdr3

    output:
    path "cdr3_sharing_pgen.tsv", emit: "shared_cdr3"
    path "sample_mapping.tsv", emit: "sample_mapping"

    script:
    """
    python - <<EOF
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    # Load data
    df = pd.read_csv("${concat_cdr3}", sep="\t")

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
        df.groupby('CDR3b')['sample_id']
        .unique()     # UNIQUE — fast & vectorized
        .apply(np.sort)  # SORT — vectorized
        .reset_index()
    )

    # Calculate counts
    grouped['total_samples'] = grouped['sample_id'].apply(len)
    grouped['samples_present'] = grouped['sample_id'].apply(
        lambda arr: ",".join(arr.astype(str))
    )

    # Drop raw list
    final_df = grouped[['CDR3b', 'total_samples', 'samples_present']]
    final_df = final_df.sort_values(by="total_samples", ascending=False)

    # Export final list
    final_df.to_csv("cdr3_sharing.tsv", sep="\t", index=False)
    EOF


    olga-compute_pgen --humanTRB -i cdr3_sharing.tsv -o pgen_sharing.tsv


    python - <<EOF
    import pandas as pd
    
    # Load TSVs for shared cdr3s and corresponding pgen values
    left_df = pd.read_csv('pgen_sharing.tsv', sep='\t', header=None, usecols=[0, 1], names=['CDR3b', 'pgen'])
    right_df = pd.read_csv('cdr3_sharing.tsv', sep='\t')

    # Drop rows where pgen == 0 and merge
    left_df = left_df[left_df['pgen'] != 0]
    merged_df = pd.merge(left_df, right_df, on='CDR3b', how='left')
    merged_df.to_csv('cdr3_sharing_pgen.tsv', sep='\t', index=False)
    EOF
    """
}

process TCRSHARING_HISTOGRAM {
    label 'process_low'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    input:
    path shared_cdr3

    output:
    path "sharing_histogram.png"

    script:
    """
    python - <<EOF
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    
    merged_df = pd.read_csv('$shared_cdr3', sep='\t')

    # Plot histogram
    sharing = merged_df['total_samples'].values

    # Create integer bin edges from 0 to max(data)
    bins = np.arange(min(sharing), max(sharing) + 2)  # +2 to include the last value as a bin edge

    plt.figure(figsize=(8, 5))
    plt.hist(sharing, bins=bins, edgecolor='black', align='left')
    plt.xticks(bins[:-1])  # whole number positions and labels
    plt.yscale('log')

    plt.xlabel('Number of Shared Samples')
    plt.ylabel('TCR Sequence Frequency (log scale)')
    plt.title('TCR Sharing Histogram')

    # Save to file
    plt.savefig("sharing_histogram.png", dpi=300, bbox_inches="tight")
    plt.close()
    EOF
    """
}

process TCRSHARING_SCATTERPLOT {
    label 'process_low'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    input:
    path shared_cdr3

    output:
    path "sharing_pgen_scatterplot.png"

    script:
    """
    python - <<EOF
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MaxNLocator

    merged_df = pd.read_csv('$shared_cdr3', sep='\t')

    # Create scatter plot with log-transform pgen
    merged_df["log10_pgen"] = np.log10(merged_df["pgen"])
    plt.figure(figsize=(8, 6))
    plt.grid(True)
    plt.scatter(merged_df["log10_pgen"], merged_df["total_samples"], c='blue', alpha=0.7)
    plt.gca().yaxis.set_major_locator(MaxNLocator(integer=True))

    plt.xlabel("log10(Probability)")
    plt.ylabel("Number of Shared Samples")
    plt.title("Scatterplot of Shared TCRs vs log10(Generation Probability)")
    plt.tight_layout()
    plt.savefig("sharing_pgen_scatterplot.png", dpi=300, bbox_inches="tight")
    plt.close()
    EOF
    """
}