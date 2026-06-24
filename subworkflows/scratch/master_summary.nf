#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MASTER_SUMMARY } from '../../modules/scratch/MASTER_SUMMARY/main.nf'

workflow MASTER_SUMMARY_SW {
    take:
    seurat_rds
    export_cells

    vdj_qc_per_sample_compact
    vdj_qc_before_after_summary
    vdj_qc_sample_sheet_resolved
    vdj_qc_clone_rank_abundance

    vdj_qc_before_after_retention_fig
    vdj_qc_pairing_bar_fig
    vdj_qc_clone_rank_abundance_fig
    vdj_qc_multiple_chains_fig

    barrier_done
    project_name

    main:
    ch_notebook = Channel.fromPath(
        "${projectDir}/modules/scratch/MASTER_SUMMARY/Master_Summary_Report.qmd",
        checkIfExists: true
    )

    MASTER_SUMMARY(
        seurat_rds,
        export_cells,
        vdj_qc_per_sample_compact,
        vdj_qc_before_after_summary,
        vdj_qc_sample_sheet_resolved,
        vdj_qc_clone_rank_abundance,
        vdj_qc_before_after_retention_fig,
        vdj_qc_pairing_bar_fig,
        vdj_qc_clone_rank_abundance_fig,
        vdj_qc_multiple_chains_fig,
        ch_notebook,
        barrier_done,
        project_name
    )

    emit:
    report_html = MASTER_SUMMARY.out.report_html
    tables      = MASTER_SUMMARY.out.tables
    figures     = MASTER_SUMMARY.out.figures
}
