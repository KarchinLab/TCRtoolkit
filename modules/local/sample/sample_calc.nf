process SAMPLE_CALC {
    tag "${sample_meta.sample}"
    label 'process_single'
    publishDir enabled: false

    input:
    tuple val(sample_meta), path(count_table)

    output:
    path "sample_stats_${sample_meta.sample}.csv"  , emit: sample_csv
    path "v_family_${sample_meta.sample}.csv"      , emit: v_family_csv
    path "d_family_${sample_meta.sample}.csv"      , emit: d_family_csv
    path "j_family_${sample_meta.sample}.csv"      , emit: j_family_csv
    val sample_meta                                  , emit: sample_meta

    script:
    """
    sample_calc.py -s '${sample_meta.sample}' -c ${count_table}
    """

    stub:
    """
    touch sample_stats.csv
    touch v_family.csv
    touch d_family.csv
    touch j_family.csv
    """
}
