/*
 * RESULTS_EXPORT
 *
 * The run's results as a few tidy, non-redundant tables in <outdir>/results/, written on
 * every run whether or not the report notebooks render (--run_reports false):
 *
 *   samples.tsv, clones.tsv.gz, clone_tests.tsv.gz, sample_overlap.tsv,
 *   cells.tsv.gz (single-cell), marker_expression.tsv (GEX), README.md
 *
 * Every report figure is a view of one of these tables (README.md maps each figure to its
 * table). Pure joins of tables the analysis modules already wrote: no Quarto, no Seurat.
 */
process RESULTS_EXPORT {
    tag "${project_name}"
    label 'process_medium'
    container "${params.container}"

    // The outputs already sit under results/, so they publish to <outdir>/results/.
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    // Module tables found by name inside in/: clonotypes.csv, *_clones_summary.csv,
    // detailed_comparisons_*.csv, sample_stats.csv, overlap_samples.tsv, olga_pgen.tsv,
    // *.vdjmatch.txt, and concatenated_cdr3_sorted.tsv as the fallback clone table.
    path tables,         stageAs: 'in/*'
    path cluster_export, stageAs: 'in_cluster_export.tsv'
    path cells,          stageAs: 'in_cells.tsv'
    path conga_cells,    stageAs: 'in_conga_cells.tsv'
    path tcri_cells,     stageAs: 'in_tcri_cells.tsv'
    path markers,        stageAs: 'in_markers.tsv'
    path samplesheet,    stageAs: 'in_samplesheet.csv'
    path script,         stageAs: 'export_results.py'
    val  project_name

    output:
    path "results/*", emit: tables

    script:
    def pcol = params.patient_col ?: 'patient'
    """
    python3 export_results.py \\
        --indir in \\
        --cluster_export in_cluster_export.tsv \\
        --cells in_cells.tsv \\
        --conga_cells in_conga_cells.tsv \\
        --tcri_cells in_tcri_cells.tsv \\
        --markers in_markers.tsv \\
        --samplesheet in_samplesheet.csv \\
        --patient_col '${pcol}' \\
        --outdir results
    """
}
