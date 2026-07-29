/*
 * Bridge: BULK_TO_EXPORT
 *
 * Builds a per-cell export table from pseudobulk clonotype data so REPERTOIRE and
 * MASTER_SUMMARY can run in the VDJ-only route (no GEX Seurat object). Each clonotype
 * is expanded into duplicate_count synthetic "cells"; annotation is constant
 * ("Unannotated") since there is no gene-expression data.
 */
process BULK_TO_EXPORT {
    tag "bulk → export_cells"
    label 'process_low'
    container "${params.container}"
    publishDir "${params.outdir}/bridge/bulk_to_export", mode: 'copy', overwrite: true

    input:
    path concat_cdr3
    path samplesheet

    output:
    path "export_cells.tsv", emit: export_cells

    script:
    def ss = (samplesheet && samplesheet.name != 'NO_FILE') ? "--samplesheet ${samplesheet}" : ""
    """
    bulk_to_export.py ${concat_cdr3} ${ss} --out export_cells.tsv
    """
}
