// process to plot sample level statistics
process SAMPLE_PLOT {
    tag "${sample_stats_csv}"
    label 'process_single'
    container "ghcr.io/karchinlab/tcrtoolkit:main"
    
    input:
    path sample_table
    path sample_stats_template
    path sample_stats_csv
    path v_family_csv

    output:
    path 'sample_stats.html'

    script:
    """
    ## copy quarto notebook to output directory
    cp $sample_stats_template sample_stats.qmd

    ## render qmd report to html
    quarto render sample_stats.qmd \
        -P project_name:$params.project_name \
        -P workflow_cmd:'$workflow.commandLine' \
        -P project_dir:$projectDir \
        -P sample_table:$sample_table \
        -P sample_stats_csv:$sample_stats_csv \
        -P v_family_csv:$v_family_csv \
        -P samplechart_x_col:${params.samplechart_x_col} \
        -P samplechart_color_col:${params.samplechart_color_col} \
        -P vgene_subject_col:${params.vgene_subject_col} \
        -P vgene_x_cols:${params.vgene_x_cols} \
        --to html
    """

    stub:
    """
    echo "1"
    """
}
