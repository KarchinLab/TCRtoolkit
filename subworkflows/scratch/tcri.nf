#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { TCRI } from '../../modules/scratch/TCRI/main.nf'

workflow TCRI_SW {
    take:
    seurat_rds
    export_cells
    project_name

    main:
    ch_notebook = Channel.fromPath(
        "${projectDir}/modules/scratch/TCRI/TCRi_Report.qmd",
        checkIfExists: true
    )

    ch_tcri = TCRI( ch_notebook, seurat_rds, export_cells, project_name )

    emit:
    report_html      = ch_tcri.report_html
    seurat_with_tcri = ch_tcri.seurat_with_tcri
    export_cells     = ch_tcri.export_cells
}
