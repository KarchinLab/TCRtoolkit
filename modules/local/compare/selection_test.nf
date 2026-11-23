/*
 * SELECTION_TEST
 *
 * Tests whether each clonotype is shared by more subjects than V(D)J generation alone
 * predicts (bin/selection_test.py). Runs on CDR3beta in every route.
 */
process SELECTION_TEST {
    label 'process_low'

    input:
    path concat_cdr3_pgen
    path samplesheet
    path anchors

    output:
    path "selection_test.tsv",        emit: results
    path "selection_summary.tsv",     emit: summary
    path "selection_calibration.tsv", emit: calibration

    script:
    """
    selection_test.py \\
        --clones ${concat_cdr3_pgen} \\
        --samplesheet ${samplesheet} \\
        --anchors ${anchors} \\
        --subject-col '${params.subject_col}' \\
        --fdr ${params.selection_fdr} \\
        --min-subjects ${params.selection_min_subjects}
    """
}
