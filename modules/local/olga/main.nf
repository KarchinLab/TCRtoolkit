process OLGA_CALCULATE{
    label 'process_single'
    publishDir enabled: false

    input:
    path cdr3_chunk

    output:
    path "pgen_${cdr3_chunk}", emit: pgen_chunk

    script:
    """
    olga-compute_pgen --humanTRB -i ${cdr3_chunk} -o olga_pgen.tsv
    cat olga_pgen.tsv \
    | awk -F'\t' -v OFS='\t' '
        BEGIN { log10_base = log(10) }
        {
            cdr3 = \$1
            pgen = \$2

            if (pgen == "NA" || pgen == "" || pgen == 0) {
                print cdr3, pgen, "NA"
            } else {
                print cdr3, pgen, log(pgen)/log10_base
            }
        }
    ' > "pgen_${cdr3_chunk}"
    """
}

process OLGA_CONCATENATE {
    label 'process_single'

    input:
    path pgen_chunks

    output:
    path "olga_pgen.tsv", emit: cdr3_pgen
    path "olga_pgen.stats.tsv", emit: cdr3_pgen_stats

    script:
    """
    awk -F'\t' '
        BEGIN {
            max_len = 0
            min_log = ""
            max_log = ""
        }

        {
            cdr3 = \$1
            log10 = \$3

            len = length(cdr3)
            if (len > max_len) max_len = len

            if (log10 != "NA" && log10 != "") {
                if (min_log == "" || log10 < min_log) min_log = log10
                if (max_log == "" || log10 > max_log) max_log = log10
            }
        }

        END {
            printf "max_cdr3_length\t%d\\n", max_len > "olga_pgen.stats.tsv"
            printf "min_log10_pgen\t%s\\n", min_log > "olga_pgen.stats.tsv"
            printf "max_log10_pgen\t%s\\n", max_log > "olga_pgen.stats.tsv"
        }
    ' ${pgen_chunks}

    printf "CDR3b\tpgen\tlog10_pgen\n" > olga_pgen.tsv
    cat ${pgen_chunks} >> olga_pgen.tsv
    """
}

process OLGA_MERGE {
    label 'process_single'

    input:
    path sorted_concat_cdr3
    path cdr3_pgen

    output:
    path "concatenated_cdr3_pgen.tsv", emit: concat_cdr3_pgen

    script:
    """
    head -n 1 ${sorted_concat_cdr3} \
        | sed 's/\$/\tpgen\tlog10_pgen/' \
        > concatenated_cdr3_pgen.tsv
    
    join -t \$'\t' \
        -1 1 -2 1 \
        -a 1 \
        -e "NA" \
        -o auto \
        <(tail -n +2 ${sorted_concat_cdr3}) \
        <(tail -n +2 ${cdr3_pgen}) \
    >> concatenated_cdr3_pgen.tsv
    """
}

process OLGA_SAMPLE_MERGE {
    tag "${sample_meta.sample}"
    label 'process_low'

    input:
    tuple val(sample_meta), path(count_table)
    path cdr3_pgen
    val olga_stats

    output:
    tuple val(sample_meta), path("${sample_meta.sample}_tcr_pgen.tsv"), emit: olga_pgen

    script:
    """
    # Extract vector of cdr3 aa, dropping null values
    python - <<EOF
import numpy as np
import pandas as pd

INDEX_FILE = "${cdr3_pgen}"
REP_FILE = "${count_table}"
OUTPUT_FILE = "${sample_meta.sample}_tcr_pgen.tsv"

CHUNKSIZE = 2_000_000


def load_index():

    print("Loading pgen index...")

    df = pd.read_csv(
        INDEX_FILE,
        sep="\t",
        usecols=["CDR3b", "pgen", "log10_pgen"],
        dtype={"CDR3b": "object"},
    )

    max_len = "${olga_stats.max_cdr3_length}"
    cdr3_dtype = f"S{max_len}"

    print(f"Max CDR3 length: {max_len}")

    cdr3_index = df["CDR3b"].values.astype(cdr3_dtype)

    pgen_index = df["pgen"].to_numpy(dtype=np.float64)

    log10_index = pd.to_numeric(
        df["log10_pgen"], errors="coerce"
    ).to_numpy(dtype=np.float32)

    del df

    print("Index loaded:", len(cdr3_index))

    return cdr3_index, pgen_index, log10_index, cdr3_dtype


def annotate_chunk(chunk, cdr3_index, pgen_index, log10_index, cdr3_dtype):

    query = chunk["junction_aa"].values.astype(cdr3_dtype)

    pos = np.searchsorted(cdr3_index, query)

    pos_safe = pos.clip(max=len(cdr3_index) - 1)

    match = (pos < len(cdr3_index)) & (cdr3_index[pos_safe] == query)

    pgen_out = np.full(len(chunk), np.nan, dtype=np.float64)
    log10_out = np.full(len(chunk), np.nan, dtype=np.float32)

    pgen_out[match] = pgen_index[pos[match]]
    log10_out[match] = log10_index[pos[match]]

    chunk["pgen"] = pgen_out
    chunk["log10_pgen"] = log10_out

    return chunk


def stream_join():

    cdr3_index, pgen_index, log10_index, cdr3_dtype = load_index()

    reader = pd.read_csv(
        REP_FILE,
        sep="\t",
        chunksize=CHUNKSIZE,
        dtype={"junction_aa": "object"},
        usecols=[
            "sequence_id",
            "junction_aa",
            "duplicate_count",
            "duplicate_frequency_percent",
        ])

    header_written = False

    for i, chunk in enumerate(reader):

        print(f"Processing chunk {i} ({len(chunk)} rows)")
        chunk = chunk[chunk['junction_aa'].notna()]

        chunk = annotate_chunk(
            chunk,
            cdr3_index,
            pgen_index,
            log10_index,
            cdr3_dtype,
        )

        chunk.to_csv(
            OUTPUT_FILE,
            sep="\t",
            index=False,
            mode="a",
            header=not header_written,
        )

        header_written = True

    print("Join complete")


if __name__ == "__main__":
    stream_join()

EOF
    """
}

process OLGA_HISTOGRAM_CALC {
    tag "${sample_meta.sample}"
    label 'process_low'

    input:
    tuple val(sample_meta), path(olga_pgen)
    val olga_stats

    output:
    tuple val(sample_meta), path("${sample_meta.sample}_histogram_data.hdf5"), emit: olga_histogram
    path "olga_ymax_value.txt", emit: olga_ymax

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

    min_val = ${olga_stats.min_log10_pgen}
    max_val = ${olga_stats.max_log10_pgen}
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