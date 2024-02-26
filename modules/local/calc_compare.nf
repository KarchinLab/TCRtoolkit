process CALC_COMPARE {
    // tag "${sample_utf8}"
    label 'process_single'

    // beforeScript 'export DOCKER_OPTS="-v $${params.data_dir}:$${params.data_dir}"'

    container "domebraccia/bulktcr:1.0-beta"

    input:
    // tuple val(sample_meta), path(count_table)
    path sample_utf8
    path meta_data

    output:
    path 'jaccard_mat.csv', emit: jaccard_mat
    path 'sorensen_mat.csv', emit: sorensen_mat
    path 'morisita_mat.csv', emit: morisita_mat

    script:
    """
    python $projectDir/bin/calc_compare.py \
        -s $sample_utf8 \
        -m $meta_data \
        -p $projectDir 
    """

}