// Generic process to render a Quarto notebook to HTML
process RENDER_NOTEBOOK {
    tag "${notebook.getBaseName()}"
    label 'process_single'

    input:
    tuple path(notebook), path(files) // path(files) just stages files in root dir
    val project_name
    val workflow_cmd

    output:
    path "${notebook.getBaseName()}.html"

    script:
    """
    ## render qmd report to html
    quarto render ${notebook} \\
        -P project_name:${project_name} \\
        -P workflow_cmd:'${workflow_cmd}' \\
        -P sample_table:${file(params.samplesheet)} \\
        --to html
    """

    stub:
    """
    touch ${notebook.getBaseName()}.html
    """
}
