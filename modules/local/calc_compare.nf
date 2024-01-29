process CALC_COMPARE {
    tag "${sample_utf8}"
    label 'process_single'

    // beforeScript 'export DOCKER_OPTS="-v $${params.data_dir}:$${params.data_dir}"'

    container "domebraccia/bulktcr:1.0-beta"

    input:
    // tuple val(sample_meta), path(count_table)
    path sample_utf8
    path meta_data

    output:
    path 'jaccard_amat.csv', emit: jaccard_amat
    path 'sorensen_amat.csv', emit: sorensen_amat
    path 'morisita_amat.csv', emit: morisita_amat

    script:
    """
    python $projectDir/bin/calc_compare.py \
        -s $sample_utf8 \
        -m $meta_data
    """

}