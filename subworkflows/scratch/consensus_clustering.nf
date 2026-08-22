#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { CONSENSUS_CLUSTERING } from '../../modules/scratch/CONSENSUS_CLUSTERING/main.nf'

workflow CONSENSUS_SW {
    take:
    seurat_rds
    export_cells
    gliph_export_cells
    tcrdist_export_cells
    giana_export_cells
    project_name

    main:
    ch_notebook = channel.fromPath(
        "${projectDir}/modules/scratch/CONSENSUS_CLUSTERING/Clonotype_Clustering_Consensus_Report.qmd",
        checkIfExists: true
    )

    CONSENSUS_CLUSTERING(
        seurat_rds,
        export_cells,
        gliph_export_cells,
        tcrdist_export_cells,
        giana_export_cells,
        ch_notebook,
        project_name
    )

    emit:
    report_html           = CONSENSUS_CLUSTERING.out.report_html
    seurat_with_consensus = CONSENSUS_CLUSTERING.out.seurat_with_consensus
    export_cells          = CONSENSUS_CLUSTERING.out.export_cells
    tables                = CONSENSUS_CLUSTERING.out.tables
    figures               = CONSENSUS_CLUSTERING.out.figures
}
