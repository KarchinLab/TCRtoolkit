#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { VDJ_QC } from '../../modules/scratch/VDJ_QC/main.nf'

workflow VDJ_QC_SW {
    take:
    ch_sample_sheet
    ch_project_name
    ch_input_annotated_object

    main:
    ch_notebook = channel.fromPath(
        "${projectDir}/modules/scratch/VDJ_QC/VDJ_QC_analysis.qmd",
        checkIfExists: true
    )

    ch_vdj_qc = VDJ_QC(
        ch_notebook,
        ch_sample_sheet,
        ch_input_annotated_object,
        ch_project_name
    )

    emit:
    report_html      = ch_vdj_qc.report_html
    contigs_after_qc = ch_vdj_qc.contigs_after_qc
    qc_tables        = ch_vdj_qc.tables
    qc_figures       = ch_vdj_qc.figures
}
