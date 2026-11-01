#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MASTER_SUMMARY } from '../../modules/scratch/MASTER_SUMMARY/main.nf'

/*
 * MASTER_SUMMARY_SW
 *
 * Runs on every route. Each `*_tables` input is the collected tables/ directory of one
 * upstream module, or an empty list on a route where that module did not run - the report
 * then explains the absence instead of rendering a blank section. `discovery_md` is the
 * Discovery Brief flattened to markdown (REPORT.out.discovery_md), shown as the report's
 * opening sections, or NO_FILE without the report notebooks.
 */
workflow MASTER_SUMMARY_SW {
    take:
    seurat_rds
    export_cells

    rollup_tables
    tcell_tables
    conga_tables
    consensus_tables
    tcri_tables
    clonal_tables       // CLONAL_DYNAMICS tables (cell states of expanding/contracting clones)
    discovery_md

    barrier_done
    project_name

    main:
    // The report plus its partials and the standalone GEX wrapper, staged together.
    ch_notebook = channel.fromPath(
        "${projectDir}/modules/scratch/MASTER_SUMMARY/*.qmd",
        checkIfExists: true
    ).collect()

    MASTER_SUMMARY(
        seurat_rds,
        export_cells,
        rollup_tables,
        tcell_tables,
        conga_tables,
        consensus_tables,
        tcri_tables,
        clonal_tables,
        discovery_md,
        ch_notebook,
        barrier_done,
        project_name
    )

    emit:
    report_html = MASTER_SUMMARY.out.report_html
    tables      = MASTER_SUMMARY.out.tables
    figures     = MASTER_SUMMARY.out.figures
    gex_html    = MASTER_SUMMARY.out.gex_html
}
