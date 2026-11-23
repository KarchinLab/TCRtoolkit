/*
 * CLONAL_PLASTICITY
 *
 * How far each clone's measured cell-state distribution moves between successive timepoints
 * (bin/clonal_plasticity.py). VDJ+GEX route only: it needs a cell state per cell.
 */
process CLONAL_PLASTICITY {
    tag "${params.project_name}"
    label 'process_medium'
    container "${params.container}"

    publishDir "${params.outdir}/clonal/plasticity", mode: 'copy', overwrite: true

    input:
    path export_cells
    path samplesheet

    output:
    path "clonal_plasticity.tsv", emit: table

    script:
    """
    clonal_plasticity.py \\
        --cells ${export_cells} \\
        --samplesheet ${samplesheet} \\
        --timepoint-order '${params.timepoint_order}' \\
        --min-cells ${params.plasticity_min_cells} \\
        --n-boot ${params.plasticity_n_boot} \\
        --n-perm ${params.plasticity_n_perm}
    """
}
