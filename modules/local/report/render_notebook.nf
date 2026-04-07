// Generic process to render a Quarto notebook to HTML
process RENDER_NOTEBOOK {
    tag "${report_name}"
    label 'process_single'

    input:
    tuple val(report_name), path(notebook), val(data_dir)
    val project_name
    val workflow_cmd

    output:
    path "${report_name}.html"

    script:
    """
    ## copy quarto notebook to working directory
    cp $notebook ${report_name}.qmd

    ## render qmd report to html
    quarto render ${report_name}.qmd \\
        -P project_name:${project_name} \\
        -P workflow_cmd:'${workflow_cmd}' \\
        -P project_dir:${data_dir} \\
        -P sample_table:${file(params.samplesheet)} \\
        --to html
    """

    stub:
    """
    touch ${report_name}.html
    """
}
