process CALC_SAMPLE {
    tag "${sample_meta[1]}"
    label 'process_single'

    container "domebraccia/bulktcr:1.0-beta"

    // publishDir "${params.output_dir}/sample_calc", mode: 'copy'

    input:
    tuple val(sample_meta), path(count_table)
    path meta_data

    output:
    path 'sample_stats.csv', emit: sample_csv
    path 'v_family.csv', emit: v_family_csv
    path 'd_family.csv', emit: d_family_csv
    path 'j_family.csv', emit: j_family_csv
    val sample_meta        , emit: sample_meta

    script:
    """
    python $projectDir/bin/calc_sample.py \
        -s '${sample_meta}' \
        -c ${count_table} \
        -m ${meta_data} 
    """

    stub:
    """
    touch sample_calc.csv
    touch v_family.csv
    touch d_family.csv
    touch j_family.csv
    """
}
