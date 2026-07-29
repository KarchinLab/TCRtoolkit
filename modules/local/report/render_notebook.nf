// Generic process to render a Quarto notebook to HTML
process RENDER_NOTEBOOK {
    tag "${notebook.getBaseName()}"
    label 'process_single'

    input:
    tuple path(notebook), path(files) // path(files) just stages files in root dir
    val project_name
    val workflow_cmd

    output:
    path "${notebook.getBaseName()}.html", emit: report_html

    script:
    // Bulk mode passes --samplesheet; single-cell mode passes --sample_sheet. Fall back so the
    // shared QC report renders in both.
    def sample_table = file(params.samplesheet ?: params.sample_sheet)
    """
    ## render qmd report to html
    quarto render ${notebook} \\
        -P project_name:${project_name} \\
        -P workflow_cmd:'${workflow_cmd}' \\
        -P sample_table:${sample_table} \\
        --to html
    """

    stub:
    """
    touch ${notebook.getBaseName()}.html
    """
}
