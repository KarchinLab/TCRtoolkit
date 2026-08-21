#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { VDJ_QC } from '../../modules/scratch/VDJ_QC/main.nf'

// A workflow-local `def as_file_list = { ... }` closure isn't callable from
// this Nextflow version's strict-syntax parser ("as_file_list is not
// defined") - a plain top-level function is, matching the same workaround
// already used for `enabled()` in workflows/tcrtoolkit_sc.nf.
def as_file_list(param_value, nofile) {
    (param_value ?: '').split(',', -1).collect { entry -> entry ? file(entry) : nofile }
}

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

    // Cirro/S3-backed datasets can't stage a whole outs/ directory - only
    // individual files (see .cirro_singlecell_*/preprocess.py, which sets
    // these three params). split(',', -1) keeps trailing empty entries so
    // list position always lines up 1:1 with sample_sheet's own row order,
    // even when a sample's file is missing; empty entries and the unset
    // (local/on-prem) case both resolve to the NO_FILE sentinel already used
    // elsewhere in this workflow.
    def nofile = file("${projectDir}/assets/NO_FILE")

    contigs_files    = as_file_list(params.contigs_files, nofile)
    clonotypes_files = as_file_list(params.clonotypes_files, nofile)
    metrics_files    = as_file_list(params.metrics_files, nofile)

    ch_vdj_qc = VDJ_QC(
        ch_notebook,
        ch_sample_sheet,
        ch_input_annotated_object,
        ch_project_name,
        contigs_files,
        clonotypes_files,
        metrics_files
    )

    emit:
    report_html      = ch_vdj_qc.report_html
    contigs_after_qc = ch_vdj_qc.contigs_after_qc
    qc_tables        = ch_vdj_qc.tables
    qc_figures       = ch_vdj_qc.figures
}
