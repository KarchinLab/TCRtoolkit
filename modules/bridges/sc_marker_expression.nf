/*
 * SC_MARKER_EXPRESSION
 *
 * Mean expression of the T-cell marker panel per cell state and clone group (expanded /
 * singleton / no TCR), from the annotated Seurat object. The one single-cell report input
 * that needs the expression matrix, computed once here so MASTER_SUMMARY renders from a
 * small table and never loads the object. Also published in results/ by RESULTS_EXPORT.
 */
process SC_MARKER_EXPRESSION {
    tag "${project_name}"
    container "${params.sc_container}"

    cpus 2
    // An RDS expands several-fold in memory; size the task from the file.
    memory {
        def gb = seurat_rds.size() / 1024**3
        gb > 40 ? 512.GB * task.attempt :
        gb > 16 ? 256.GB * task.attempt :
        gb >  6 ? 128.GB * task.attempt :
        gb >  2 ?  64.GB * task.attempt :
                   16.GB * task.attempt
    }

    publishDir "${params.outdir}/sc/marker_expression", mode: 'copy', overwrite: true

    input:
    path seurat_rds,   stageAs: 'in_seurat.rds'
    path export_cells, stageAs: 'in_export_cells.tsv'
    // Staged rather than called from bin/ so edits to it change the task hash.
    path script,       stageAs: 'sc_marker_expression.R'
    val  project_name

    output:
    path "marker_expression.tsv", emit: table

    script:
    """
    Rscript sc_marker_expression.R \\
        --seurat_rds in_seurat.rds \\
        --export_cells in_export_cells.tsv \\
        --out marker_expression.tsv
    """
}
