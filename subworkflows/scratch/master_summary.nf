#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MASTER_SUMMARY } from '../../modules/scratch/MASTER_SUMMARY/main.nf'

/*
 * MASTER_SUMMARY_SW
 *
 * Runs on both single-cell routes. Each `*_tables` input is the collected tables/
 * directory of one upstream module, or an empty list on a route where that module did
 * not run - the report then explains the absence instead of rendering a blank section.
 */
workflow MASTER_SUMMARY_SW {
    take:
    seurat_rds
    export_cells

    vdj_qc_tables
    pseudobulk_tables
    sample_tables
    compare_tables
    rollup_tables
    repertoire_tables
    tcell_tables
    conga_tables
    consensus_tables
    tcri_tables
    mergedvdj_tables

    barrier_done
    project_name

    main:
    ch_notebook = channel.fromPath(
        "${projectDir}/modules/scratch/MASTER_SUMMARY/Master_Summary_Report.qmd",
        checkIfExists: true
    )

    MASTER_SUMMARY(
        seurat_rds,
        export_cells,
        vdj_qc_tables,
        pseudobulk_tables,
        sample_tables,
        compare_tables,
        rollup_tables,
        repertoire_tables,
        tcell_tables,
        conga_tables,
        consensus_tables,
        tcri_tables,
        mergedvdj_tables,
        ch_notebook,
        barrier_done,
        project_name
    )

    emit:
    report_html = MASTER_SUMMARY.out.report_html
    tables      = MASTER_SUMMARY.out.tables
    figures     = MASTER_SUMMARY.out.figures
}
