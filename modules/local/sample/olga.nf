process OLGA_PGEN_CALC {
    tag "${sample_meta.sample}"
    label 'process_low'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    input:
    tuple val(sample_meta), path(count_table)

    output:
    tuple val(sample_meta), path("${sample_meta.sample}_tcr_pgen.tsv"), emit: "olga_pgen"
    path "olga_xmin_value.txt", emit: 'olga_xmin'
    path "olga_xmax_value.txt", emit: 'olga_xmax'

    script:
    """
    # Extract vector of cdr3 aa, dropping null values
    python - <<EOF
    import pandas as pd

    df = pd.read_csv("${count_table}", sep="\t")
    df = df.dropna(subset=["junction_aa"])
    df = df["junction_aa"]
    df.to_csv("output.tsv", sep="\t", index=False, header=False)
    EOF


    olga-compute_pgen --humanTRB -i output.tsv -o "${sample_meta.sample}_pgen.tsv"


    python - <<EOF
    import numpy as np
    import pandas as pd

    # Merge count and probability generation tables 
    df1 = pd.read_csv("${count_table}", sep="\t")
    df1 = df1.dropna(subset=["junction_aa"])
    df2 = pd.read_csv('${sample_meta.sample}_pgen.tsv', sep='\t', header=None, usecols=[0, 1], names=['junction_aa', 'pgen'])
    merged_df = pd.merge(df1, df2, on='junction_aa', how='left')
    merged_df.to_csv("${sample_meta.sample}_tcr_pgen.tsv", sep="\t", index=False)

    merged_df = merged_df[merged_df['pgen'] != 0]
    log_probs = np.log10(merged_df['pgen'])

    left_bound = np.floor(np.min(log_probs))
    right_bound = np.ceil(np.max(log_probs))

    with open(f"olga_xmin_value.txt", "w") as f:
        f.write(str(left_bound))

    with open(f"olga_xmax_value.txt", "w") as f:
        f.write(str(right_bound))
    EOF
    """
}

process OLGA_HISTOGRAM_CALC {
    tag "${sample_meta.sample}"
    label 'process_low'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    input:
    tuple val(sample_meta), path(olga_pgen)
    val olga_global_xmin
    val olga_global_xmax

    output:
    tuple val(sample_meta), path("${sample_meta.sample}_histogram_data.hdf5"), emit: "olga_histogram"
    path "olga_ymax_value.txt", emit: 'olga_ymax'

    script:
    """
    python - <<EOF
    import h5py
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    merged_df = pd.read_csv("${olga_pgen}", sep="\t")

    # Drop rows where pgen is 0
    merged_df = merged_df[merged_df['pgen'] != 0]
    log_probs = np.log10(merged_df['pgen'])
    cdr3_counts = merged_df['duplicate_count']

    min_val = ${olga_global_xmin}
    max_val = ${olga_global_xmax}
    n_bins = 30

    bin_edges = np.linspace(min_val, max_val, n_bins + 1)

    # Manually compute weighted histogram
    counts = np.zeros(len(bin_edges) - 1)
    for log_p, weight in zip(log_probs, cdr3_counts):
        bin_idx = np.searchsorted(bin_edges, log_p, side='right') - 1
        if 0 <= bin_idx < len(counts):
            counts[bin_idx] += weight

    # Save histogram data
    with h5py.File(f"${sample_meta.sample}_histogram_data.hdf5", "w") as f:
        f.create_dataset("counts", data=counts)
        f.create_dataset("bin_edges", data=bin_edges)

    # Save max count value (rounded up) for y-axis standardization
    normalized_counts = counts / counts.sum()
    ymax = np.ceil(np.max(normalized_counts) * 100) / 100
    with open("olga_ymax_value.txt", "w") as f:
        f.write(f"{ymax}\\n")
    EOF
    """
}


process OLGA_HISTOGRAM_PLOT {
    tag "${sample_meta.sample}"
    label 'process_low'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    input:
    tuple val(sample_meta), path(olga_histogram)
    val olga_global_ymax

    output:
    path "${sample_meta.sample}_tcr_pgen_histogram.png"

    script:
    """
    python - <<EOF
    import h5py
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    # Load histogram data
    with h5py.File(f"${olga_histogram}", "r") as f:
        counts = f["counts"][:]
        bin_edges = f["bin_edges"][:]
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_widths = np.diff(bin_edges)
    normalized_counts = counts / counts.sum()

    # Filter: exclude zero-count bins
    nonzero_mask = counts > 0
    filtered_centers = bin_centers[nonzero_mask]
    filtered_heights = normalized_counts[nonzero_mask]
    filtered_widths = bin_widths[nonzero_mask]

    # Plot histogram
    plt.figure(figsize=(8, 5))
    plt.bar(
        filtered_centers,
        filtered_heights,
        width=filtered_widths,
        edgecolor='black',
        align='center'
    )

    plt.xlim(bin_edges[0] - 0.5, bin_edges[-1] + 0.5)
    standardize_y = True # Could set to a future param if desired
    if standardize_y:
        plt.ylim(top=${olga_global_ymax})

    plt.xlabel('log_10 Generation Probability')
    plt.ylabel('Probability Density')
    plt.title(f'${sample_meta.sample} TCR Generation Probability Histogram')

    # Save to file
    plt.savefig("${sample_meta.sample}_tcr_pgen_histogram.png", dpi=300, bbox_inches="tight")
    plt.close()
    EOF
    """
}

process OLGA_WRITE_MAX {
    label 'process_single'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    input:
    val olga_global_xmin
    val olga_global_xmax
    val olga_global_ymax

    output:
    path "olga_xmin_value.txt"
    path "olga_xmax_value.txt"
    path "olga_ymax_value.txt"

    script:
    """
    echo ${olga_global_xmin} > olga_xmin_value.txt
    echo ${olga_global_xmax} > olga_xmax_value.txt
    echo ${olga_global_ymax} > olga_ymax_value.txt
    """
}