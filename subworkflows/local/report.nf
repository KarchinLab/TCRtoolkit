/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RENDER_NOTEBOOK } from '../../modules/local/report/render_notebook'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Each main notebook ends with {{< include ./_sc_<name>.qmd >}}: the single-cell-only
// section appended to it (CoNGA, TCRi, phenotype, VDJ QC, cell-level clustering).
def scSection(notebookName) {
    [
        template_qc              : '_sc_qc.qmd',
        template_discovery_brief : '_sc_discovery_brief.qmd',
        template_details_sample  : '_sc_details_sample.qmd',
        template_details_compare : '_sc_details_compare.qmd',
        template_details_patient : '_sc_details_patient.qmd',
    ][notebookName]
}

// Single-channel listener: receives tuples of
// (notebook_template, files_to_stage, staged_layout)
// and renders each notebook to HTML.
workflow REPORT {

    take:
    ch_reports  // channel of tuples: tuple(path notebook_template, path files_to_stage, staged_layout)
    samplesheet // value channel: samplesheet the notebooks read as sample_table
    sc_bundle   // value channel: SC_REPORT_BUNDLE's sc_bundle/ directory, or NO_FILE (bulk)

    main:
    def notebooks_dir = file(params.template_qc).parent

    // Quarto includes are static, so the SC section is always included: with SC data the
    // real partial is staged alongside the sc_bundle/ tables it reads; without, an empty
    // stub is staged under the partial's name and the notebook renders as before.
    ch_with_sc = ch_reports
        .combine(sc_bundle)
        .map { notebook, files, layout, bundle ->
            def section = scSection(notebook.baseName)
            if (!section) return [notebook, files, layout]
            if (bundle.name != 'NO_FILE') {
                def sc_files = [file("${notebooks_dir}/${section}"),
                                file("${notebooks_dir}/_sc_common.qmd"),
                                bundle]
                [notebook, files + sc_files, layout]
            } else {
                [notebook, files + [file("${notebooks_dir}/_sc_stub.qmd")],
                 layout + [[section, '_sc_stub.qmd']]]
            }
        }

    RENDER_NOTEBOOK(
        ch_with_sc,
        params.project_name,
        workflow.commandLine,
        samplesheet
    )

    emit:
        report_html = RENDER_NOTEBOOK.out.report_html
        figures     = RENDER_NOTEBOOK.out.figures
}
