process PLOT_COMPARE {
    tag "${sample_stats_csv}"
    label 'plot_compare'

    container "domebraccia/bulktcr:1.0-beta"

    publishDir "${params.output_dir}/plot_compare", mode: 'copy'
    
    input:
    // path sample_table
    // path sample_stats_csv

    output:
    path 'compare_stats.html'

    script:    
    """
    ## copy quarto notebook to output directory
    cp $projectDir/notebooks/compare_stats_template.qmd compare_stats.qmd

    ## render qmd report to html
    quarto render compare_stats.qmd \
        -P project_name:$params.project_name \
        -P workflow_cmd:'$workflow.commandLine' \
        -P project_dir:$projectDir \
        -P sample_table:$sample_table \
        -P sample_stats_csv:$sample_stats_csv \
        --to html
    """

    stub:
    """
    touch compare_stats.qmd
    """
    
    }