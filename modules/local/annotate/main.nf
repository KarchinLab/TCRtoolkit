process ANNOTATE_PROCESS {
    tag "${sample_meta.sample}"
    label 'process_low'
    publishDir enabled: false

    input:
    tuple val(sample_meta), path(count_table)

    output:
    tuple val(sample_meta), path("${sample_meta.sample}_cdr3.tsv"), emit: "process"
    tuple val(sample_meta), path("${sample_meta.sample}_pre_filter_stats.csv"), emit: "pre_filter_stats"

    script:
    """
    python - <<EOF
import pandas as pd

USECOLS = [
    "junction_aa", "v_call", "d_call", "j_call",
    "duplicate_count", "junction_aa_length", "duplicate_frequency_percent",
    "productive", "sequence", "sequence_id", "junction",
]

df = pd.read_csv(
    "${count_table}",
    sep='\t',
    usecols=USECOLS,
    dtype={
        "junction_aa": "string",
        "v_call": "string",
        "d_call": "string",
        "j_call": "string",
        "sequence": "string",
        "sequence_id": "string",
        "junction": "string",
        "duplicate_count": "int",
        "productive": "boolean"
    })

total_clones = len(df)
num_prod = int(df['productive'].sum())
num_nonprod = total_clones - num_prod

pd.DataFrame({
    "sample": ["${sample_meta.sample}"],
    "total_clones": [total_clones],
    "productive_clones": [num_prod],
    "nonproductive_clones": [num_nonprod],
}).to_csv("${sample_meta.sample}_pre_filter_stats.csv", index=False)

df = df[df.junction_aa.notna()]
df = df[df['productive']==True]

OUTPUT_COLS = [
    "junction_aa", "v_call", "d_call", "j_call",
    "duplicate_count", "junction_aa_length", "duplicate_frequency_percent",
    "sequence", "sequence_id", "junction",
]
df["sample"] = "${sample_meta.sample}"
df[OUTPUT_COLS + ["sample"]].to_csv("${sample_meta.sample}_cdr3.tsv", sep="\t", index=False)

EOF
    """
}

process ANNOTATE_CONCATENATE {
    label 'process_low'

    input:
    path samplesheet_utf8
    path all_sample_files

    output:
    path "concatenated_cdr3.tsv", emit: concat_cdr3

    script:
    """
    # Concatenate input Adaptive files and process metadata
    # Note: 'all_sample_files' is used as an implicit dependency to control scheduling.
    : $all_sample_files
    compare_concatenate.py "${samplesheet_utf8}"
    """
}

process ANNOTATE_SORT_CDR3 {
    label 'process_medium'

    input:
    path concat_cdr3

    output:
    path 'concatenated_cdr3_sorted.tsv', emit: concat_cdr3_sorted

    script:
    """
    head -n 1 ${concat_cdr3} > concatenated_cdr3_sorted.tsv

    tail -n +2 ${concat_cdr3} \
        | LC_ALL=C sort \
            -t \$'\t' \
            -k1,1 -k2,5 \
            --parallel=${task.cpus} \
            -S 50% \
        >> concatenated_cdr3_sorted.tsv
    """
}

process ANNOTATE_DEDUPLICATE_CDR3_TRBV {
    label 'process_low'

    input:
    path concat_cdr3

    output:
    path 'unique_cdr3_trbv.tsv', emit: unique_cdr3_trbv
    path 'unique_cdr3_trbv_with_vcall.tsv', emit: unique_cdr3_trbv_with_vcall

    script:
    """
    tail -n +2 ${concat_cdr3} \
        | awk -F'\t' '{print toupper(\$1) "\t" toupper(\$2)}' \
        | LC_ALL=C sort -u \
        > unique_cdr3_trbv.tsv

    # additional file with blank TRBV calls removed for GIANA
    awk -F'\t' 'NF>=2 && \$2 ~ /^TRBV/' unique_cdr3_trbv.tsv > unique_cdr3_trbv_with_vcall.tsv
    """
}

process ANNOTATE_DEDUPLICATE_CDR3 {
    label 'process_single'

    input:
    path unique_cdr3_trbv

    output:
    path 'unique_cdr3.txt', emit: unique_cdr3

    script:
    """
    cut -f1 ${unique_cdr3_trbv} \
        | LC_ALL=C sort -u \
        > unique_cdr3.txt
    """
}