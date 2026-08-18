#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { CONGA } from '../../modules/scratch/CONGA/main.nf'

workflow CONGA_SW {
    take:
    seurat_rds
    export_cells
    project_name

    main:
    ch_notebook = channel.fromPath(
        "${projectDir}/modules/scratch/CONGA/CoNGA_Report.qmd",
        checkIfExists: true
    )

    CONGA( seurat_rds, export_cells, ch_notebook, project_name )

    emit:
    report_html       = CONGA.out.report_html
    seurat_with_conga = CONGA.out.seurat_with_conga
    export_cells      = CONGA.out.export_cells
    data              = CONGA.out.data
    tables            = CONGA.out.tables
    figures           = CONGA.out.figures
}
