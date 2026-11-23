// Cross-study comparison: the run's clonotypes found in external reference clonotype sets
// (--reference_samplesheet). Exact matching on --reference_match_key; see
// bin/reference_match.py for the accepted reference formats.
process REFERENCE_MATCH {
    label 'process_low'

    input:
    path concat_cdr3_sorted
    path samplesheet
    // The reference samplesheet and every table it lists, staged side by side so its
    // relative `file` entries resolve.
    path reference_samplesheet, stageAs: 'reference/*'
    path reference_files,       stageAs: 'reference/*'

    output:
    path "reference_matches.tsv",         emit: matches
    path "reference_overlap_summary.tsv", emit: summary
    path "per_sample/*",                  emit: per_sample, optional: true

    script:
    """
    reference_match.py \\
        --concat ${concat_cdr3_sorted} \\
        --samplesheet ${samplesheet} \\
        --reference ${reference_samplesheet} \\
        --subject-col '${params.subject_col}' \\
        --match-key '${params.reference_match_key}'
    """

    stub:
    """
    touch reference_matches.tsv reference_overlap_summary.tsv
    """
}
