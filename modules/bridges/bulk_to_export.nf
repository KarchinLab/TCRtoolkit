/*
 * Bridge: BULK_TO_EXPORT
 *
 * Builds a per-cell export table from clonotype data so CLUSTER_TO_SC, CONSENSUS
 * and MASTER_SUMMARY can run without a GEX Seurat object. `level` is 'cell'
 * (one row per cell; single-cell VDJ-only route) or 'clonotype' (one row per sample x
 * clonotype; bulk, where counts are reads/templates). Annotation is constant
 * ("Unannotated") since there is no gene-expression data.
 */
process BULK_TO_EXPORT {
    tag "bulk → export_cells (${level})"
    label 'process_low'
    container "${params.container}"
    publishDir "${params.outdir}/bridge/bulk_to_export", mode: 'copy', overwrite: true

    input:
    path concat_cdr3
    path samplesheet
    val  level

    output:
    path "export_cells.tsv", emit: export_cells

    script:
    def ss = (samplesheet && samplesheet.name != 'NO_FILE') ? "--samplesheet ${samplesheet}" : ""
    """
    bulk_to_export.py ${concat_cdr3} ${ss} --level ${level} --out export_cells.tsv
    """
}
