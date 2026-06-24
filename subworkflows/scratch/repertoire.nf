#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { REPERTOIRE } from '../../modules/scratch/REPERTOIRE/main.nf'

workflow REPERTOIRE_SW {
    take:
    seurat_rds
    export_cells
    project_name

    main:
    ch_notebook = Channel.fromPath(
        "${projectDir}/modules/scratch/REPERTOIRE/Repertoire_Report.qmd",
        checkIfExists: true
    )

    REPERTOIRE( seurat_rds, export_cells, ch_notebook, project_name )

    emit:
    report_html = REPERTOIRE.out.report_html
    tables      = REPERTOIRE.out.tables
    figures     = REPERTOIRE.out.figures
}
