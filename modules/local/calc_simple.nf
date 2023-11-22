process SIMPLE_CALC {
    tag "${sample_meta[1]}"
    label 'process_single'

    container "domebraccia/bulktcr:0.2"

    publishDir "${params.output_dir}/simple_calc", mode: 'copy'

    input:
    tuple val(sample_meta), path(count_table)

    output:
    path 'simple_calc.csv', emit: simple_csv
    path 'gene_usage_*.pkl', emit: gene_usage_pkl

    script:
    """
    python $projectDir/bin/calc_simple.py \
        -m '${sample_meta}' \
        -c ${count_table} 
    """

    stub:
    """
    touch simple_calc.csv
    """
}