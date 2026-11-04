process CLONAL_DYNAMICS {
    label 'process_low'

    input:
    path concat_cdr3_sorted
    path samplesheet

    output:
    path "clonal_dynamics/clonotypes.csv",                        emit: clonotypes
    path "clonal_dynamics/detailed_comparisons_expansion.csv",    emit: detailed_expansion
    path "clonal_dynamics/detailed_comparisons_contraction.csv",  emit: detailed_contraction
    path "clonal_dynamics/expanded_clones_summary.csv",           emit: expanded_summary
    path "clonal_dynamics/contracting_clones_summary.csv",        emit: contracting_summary
    path "clonal_dynamics/expansion_counts.csv",                  emit: expansion_counts
    path "clonal_dynamics/contraction_counts.csv",                emit: contraction_counts

    script:
    """
    clonal_dynamics.py \\
        --concat_cdr3 ${concat_cdr3_sorted} \\
        --sample_table ${samplesheet} \\
        --subject_col '${params.subject_col}' \\
        --origin_col 'origin' \\
        --timepoint_col '${params.timepoint_col}' \\
        --timepoint_order_col '${params.timepoint_order_col}' \\
        --timepoint_order '${params.timepoint_order}' \\
        --alias_col '${params.alias_col}' \\
        --outdir clonal_dynamics
    """
}
