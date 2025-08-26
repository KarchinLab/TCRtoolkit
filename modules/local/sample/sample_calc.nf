process SAMPLE_CALC {
    tag "${sample_meta.sample}"
    label 'process_single'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    input:
    tuple val(sample_meta), path(count_table)

    output:
    path 'sample_stats.csv'     , emit: sample_csv
    path 'v_family.csv'         , emit: v_family_csv
    path 'd_family.csv'         , emit: d_family_csv
    path 'j_family.csv'         , emit: j_family_csv
    val sample_meta             , emit: sample_meta

    script:
    """
    echo '' > sample_stats.csv
    sample_calc.py -s '${sample_meta}' -c ${count_table}
    """

    stub:
    """
    touch sample_stats.csv
    touch v_family.csv
    touch d_family.csv
    touch j_family.csv
    """
}
