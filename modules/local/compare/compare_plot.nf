process COMPARE_PLOT {
    label 'process_single'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    input:
    path sample_utf8
    path jaccard_mat
    path sorensen_mat
    path morisita_mat
    path compare_stats_template
    val project_name
    path all_sample_files

    output:
    path 'compare_stats.html'

    script:    
    """
    ## copy quarto notebook to output directory
    cp $compare_stats_template compare_stats.qmd

    ## render qmd report to html
    # export QUARTO_DENO_EXTRA_OPTIONS=--v8-flags=--max-old-space-size=8192
    quarto render compare_stats.qmd \
        -P project_name:$project_name \
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
