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

// Single-channel listener: receives tuples of (report_name, notebook_path, data_dir)
// and renders each notebook to HTML.
workflow REPORT {

    take:
    ch_reports // channel of tuples: tuple(val report_name, path notebook, val data_dir)

    main:
    RENDER_NOTEBOOK(
        ch_reports,
        params.project_name,
        workflow.commandLine
    )

    emit:
    reports = RENDER_NOTEBOOK.out
}
