/*
 * Bridge 1: SC_TO_BULK
 *
 * Converts SCRATCH-TCR's per-cell export_cells.tsv (from T-cell integration)
 * into per-sample AIRR-format TSVs that TCRtoolkit's ANNOTATE step can consume.
 *
 * Logic:
 *   - Extracts beta chain CDR3 (junction_aa), TRBV (v_call), TRBJ (j_call) from CTaa/CTgene
 *   - Groups cells by sample, counts cells per clonotype → duplicate_count
 *   - Writes one AIRR TSV per sample to bulk_samples/
 *   - Writes a synthetic samplesheet.csv for downstream TCRtoolkit modules
 */

process SC_TO_BULK {
    tag "SC → Bulk (Bridge 1)"
    label 'process_low'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    publishDir "${params.outdir}/bridge/sc_to_bulk", mode: 'copy', overwrite: true

    input:
    path  export_cells
    val   sample_col
    val   meta_cols

    output:
    path "bulk_samples/*.tsv",       emit: bulk_tsv_files
    path "synthetic_samplesheet.csv", emit: samplesheet

    script:
    """
    sc_to_bulk.py \\
        "${export_cells}" \\
        "${sample_col}" \\
        "${meta_cols}"
    """
}
