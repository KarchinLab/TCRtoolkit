#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { TCRDIST3 } from '../../modules/scratch/TCRDIST3/main.nf'

workflow TCRDIST3_SW {
    take:
    seurat_rds
    export_cells
    project_name

    main:
    ch_notebook = Channel.fromPath(
        "${projectDir}/modules/scratch/TCRDIST3/TCRdist3_Report.qmd",
        checkIfExists: true
    )

    TCRDIST3( seurat_rds, export_cells, ch_notebook, project_name )

    emit:
    report_html          = TCRDIST3.out.report_html
    seurat_with_tcrdist3 = TCRDIST3.out.seurat_with_tcrdist3
    export_cells         = TCRDIST3.out.export_cells
    tables               = TCRDIST3.out.tables
    figures              = TCRDIST3.out.figures
}
