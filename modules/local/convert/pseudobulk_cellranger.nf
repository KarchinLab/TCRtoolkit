process PSEUDOBULK_CELLRANGER {
    tag "${sample_meta.sample}"
    label 'process_low'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    input:
    tuple val(sample_meta), path(count_table)
    path airr_schema

    output:
    tuple val(sample_meta), path("${sample_meta.sample}_pseudobulk.tsv") , emit: "cellranger_pseudobulk"

    script:
    """
    pseudobulk.py ${count_table} ${sample_meta.sample} ${airr_schema}
    """
}