#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { TCELL_INTEGRATION } from '../../modules/scratch/TCELL_INTEGRATION/main.nf'

workflow TCELL_INTEGRATION_SW {
    take:
    contigs_after_qc
    annotated_object
    project_name

    main:
    ch_notebook = Channel.fromPath(
        "${projectDir}/modules/scratch/TCELL_INTEGRATION/TCell_Integration_Report.qmd",
        checkIfExists: true
    )

    ch_tcell = TCELL_INTEGRATION(
        contigs_after_qc,
        annotated_object,
        ch_notebook,
        project_name
    )

    emit:
    report_html            = ch_tcell.report_html
    seurat_tcells_with_tcr = ch_tcell.seurat_tcells_with_tcr
    export_cells           = ch_tcell.export_cells
    tables                 = ch_tcell.tables
    figures                = ch_tcell.figures
}
