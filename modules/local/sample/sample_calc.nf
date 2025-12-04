process SAMPLE_CALC {
    tag "${sample_meta.sample}"
    label 'process_single'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    input:
    tuple val(sample_meta), path(count_table)

    output:
    path "stats/sample_stats_${sample_meta.sample}.csv"  , emit: sample_csv
    path "vdj/v_family_${sample_meta.sample}.csv"      , emit: v_family_csv
    path "vdj/d_family_${sample_meta.sample}.csv"      , emit: d_family_csv
    path "vdj/j_family_${sample_meta.sample}.csv"      , emit: j_family_csv
    val sample_meta                                  , emit: sample_meta

    script:
    def meta_json = groovy.json.JsonOutput.toJson(sample_meta)

    """
    mkdir -p stats
    mkdir -p vdj
    
    sample_calc.py -s '${meta_json}' -c ${count_table}
    """

    stub:
    """
    touch sample_stats.csv
    touch v_family.csv
    touch d_family.csv
    touch j_family.csv
    """
}
