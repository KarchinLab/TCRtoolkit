process CONVERT_ADAPTIVE {
    tag "${sample_meta.sample}"
    label 'process_low'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    input:
    tuple val(sample_meta), path(count_table)
    path airr_schema
    path imgt_lookup

    output:
    tuple val(sample_meta), path("${sample_meta.sample}_airr.tsv") , emit: "adaptive_convert"

    script:
    """
    convert.py ${count_table} ${sample_meta.sample} ${imgt_lookup} ${airr_schema}
    """
}