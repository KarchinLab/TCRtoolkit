process CONVERGENCE {
    tag "${sample_meta.sample}"
    label 'process_low'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    input:
    tuple val(sample_meta), path(count_table)

    output:
    path "${count_table.baseName}_tcr_convergence.tsv", emit: "convergence_output"
    path "${count_table.baseName}_tcr_convergence_histogram.png"

    script:
    """
    # Extract vector of cdr3 aa, dropping null values
    python - <<EOF
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    # Load your TCR data (make sure the file has 'cdr3_aa' and 'cdr3_nt' columns)
    df = pd.read_csv("${count_table}", sep="\t", usecols=["junction_aa", "sequence"])
    df = df.dropna(subset=["junction_aa"])

    # Group by amino acid sequence and count unique nucleotide sequences (convergence)
    convergence_df = (
        df.groupby("junction_aa")["sequence"]
        .nunique()
        .reset_index(name="convergence")
    )

    # Sort by convergence count, descending
    convergence_df = convergence_df.sort_values(by="convergence", ascending=False)

    # Export
    convergence_df.to_csv("${count_table.baseName}_tcr_convergence.tsv", sep="\t", index=False)
    
    # Plot histogram
    convergence = convergence_df['convergence'].values
    average_convergence = convergence_df["convergence"].mean()

    # Create integer bin edges from 0 to max(data)
    bins = np.arange(min(convergence), max(convergence) + 2)  # +2 to include the last value as a bin edge

    plt.figure(figsize=(8, 5))
    plt.hist(convergence, bins=bins, edgecolor='black', align='left')
    plt.xticks(bins[:-1])  # whole number positions and labels
    plt.yscale('log')

    plt.xlabel('TCR Convergence Number')
    plt.ylabel('TCR Convergence Frequency (log scale)')
    plt.title(f'${count_table.baseName} TCR Convergence Histogram, Average: {average_convergence:.2f}')

    # Save to file
    plt.savefig("${count_table.baseName}_tcr_convergence_histogram.png", dpi=300, bbox_inches="tight")
    plt.close()

    EOF
    """
}
