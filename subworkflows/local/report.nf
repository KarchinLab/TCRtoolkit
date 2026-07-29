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

// Single-channel listener: receives tuples of
// (notebook_template, files_to_stage)
// and renders each notebook to HTML.
workflow REPORT {

    take:
    ch_reports // channel of tuples: tuple(path notebook_template, path files_to_stage)
    main:
    RENDER_NOTEBOOK(
        ch_reports,
        params.project_name,
        workflow.commandLine
    )

    emit:
        rendered_html = RENDER_NOTEBOOK.out
}
