/*
 * Bridge 3: VDJ_TO_BULK
 *
 * Converts VDJ_QC's contigs_after_qc.tsv (per-contig rows, one per chain per cell)
 * into per-sample AIRR-format TSVs for TCRtoolkit.
 *
 * Used in VDJ-only mode when no GEX Seurat object is available and
 * TCELL_INTEGRATION is skipped.
 */

process VDJ_TO_BULK {
    tag "VDJ → Bulk (Bridge 3)"
    label 'process_low'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    publishDir "${params.outdir}/bridge/vdj_to_bulk", mode: 'copy', overwrite: true

    input:
    path  contigs_after_qc
    val   sample_col

    output:
    path "bulk_samples/*.tsv",        emit: bulk_tsv_files
    path "synthetic_samplesheet.csv", emit: samplesheet

    script:
    """
    vdj_to_bulk.py \\
        "${contigs_after_qc}" \\
        "${sample_col}"
    """
}
