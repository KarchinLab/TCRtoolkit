process PLOT_COMPARE {
    // tag "${jaccard_mat}"
    label 'plot_compare'

    container "domebraccia/bulktcr:1.0-beta"

    publishDir "${params.output_dir}/plot_compare", mode: 'copy'
    
    input:
    path sample_utf8
    path jaccard_mat
    path sorensen_mat
    path morisita_mat

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
        -P jaccard_mat:$jaccard_mat \
        -P sorensen_mat:$sorensen_mat \
        -P morisita_mat:$morisita_mat \
        -P sample_utf8:$sample_utf8 \
        --to html
    """

    stub:
    """
    touch compare_stats.qmd
    """
    
    }