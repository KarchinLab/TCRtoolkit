/*
 * REPERTOIRE_METRICS
 *
 * The single place repertoire-level metrics are computed, for both modalities: per-sample
 * diversity (merged into sample_stats.csv), pooled diversity per patient / condition /
 * timepoint / origin, clone-size bins, pairwise sample and patient overlap, and clone
 * sharing. The report notebooks and the Master Summary read these tables rather than
 * recomputing them. Bulk counts TRB clonotypes; single-cell counts paired clones in cells.
 */
process REPERTOIRE_METRICS {
    label 'process_low'
    container "${params.container}"

    publishDir "${params.outdir}/sample", mode: 'copy', overwrite: true, pattern: 'sample_stats.csv'
    publishDir "${params.outdir}/repertoire", mode: 'copy', overwrite: true, pattern: '*.tsv'

    input:
    path sample_stats, stageAs: 'in_sample_stats.csv'
    path clones,       stageAs: 'in_clones.tsv'
    path samplesheet,  stageAs: 'in_samplesheet.csv'

    output:
    path "sample_stats.csv", emit: sample_stats
    path "*.tsv",            emit: tables

    script:
    """
    repertoire_metrics.py \\
        --clones in_clones.tsv \\
        --sample-stats in_sample_stats.csv \\
        --samplesheet in_samplesheet.csv \\
        --bin-breaks '${params.clone_bin_breaks}' \\
        --bin-labels '${params.clone_bin_labels}' \\
        --outdir .
    """
}
