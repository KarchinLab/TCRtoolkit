// Clonotype sharing between the compartments of each patient (tumor, blood, CSF, ...):
// pairwise overlap metrics and the per-compartment clone frequencies behind the Master
// Summary's flow diagrams. Always writes its tables; when the samplesheet has no
// compartment column, or no patient was sampled in two compartments, they are headers only
// and compartment_overlap_status.txt says why (the report then skips the section).
process COMPARTMENT_OVERLAP {
    label 'process_low'

    input:
    path concat_cdr3_sorted
    path samplesheet

    output:
    path "compartment_overlap_summary.tsv", emit: summary
    path "compartment_flow.tsv",            emit: flow
    path "compartment_shared_clones.tsv",   emit: shared_clones
    path "compartment_overlap_status.txt",  emit: status
    path "compartment_samples.tsv",         emit: samples

    script:
    """
    compartment_overlap.py \\
        --concat ${concat_cdr3_sorted} \\
        --samplesheet ${samplesheet} \\
        --compartment-col '${params.compartment_col}' \\
        --subject-col '${params.subject_col}' \\
        --time-col '${params.cohort_time_col}' \\
        --fallback-time-col '${params.timepoint_col}' \\
        --match-key '${params.compartment_match_key}' \\
        --top-n ${params.compartment_top_n}
    """

    stub:
    """
    touch compartment_overlap_summary.tsv compartment_flow.tsv compartment_shared_clones.tsv compartment_samples.tsv
    echo "skipped: stub run" > compartment_overlap_status.txt
    """
}
