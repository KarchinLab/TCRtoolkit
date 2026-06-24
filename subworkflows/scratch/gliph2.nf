#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { GLIPH2 } from '../../modules/scratch/GLIPH2/main.nf'

workflow GLIPH2_SW {
    take:
    seurat_rds
    export_cells
    project_name

    main:
    ch_notebook = Channel.fromPath(
        "${projectDir}/modules/scratch/GLIPH2/GLIPH2_Report.qmd",
        checkIfExists: true
    )
    ch_ref = Channel.fromPath(params.gliph_reference_bundle, checkIfExists: true)

    GLIPH2( seurat_rds, export_cells, ch_notebook, ch_ref, project_name )

    emit:
    report_html        = GLIPH2.out.report_html
    seurat_with_gliph2 = GLIPH2.out.seurat_with_gliph2
    export_cells       = GLIPH2.out.export_cells
    tables             = GLIPH2.out.tables
    figures            = GLIPH2.out.figures
}
