process PSEUDOBULK_PHENOTYPE_CELLRANGER {
    tag "${sample_meta.sample}"
    label 'process_low'

    input:
    tuple val(sample_meta), path(count_table)
    path airr_schema
    path sobject_gex

    output:
    tuple val(sample_meta), path("${sample_meta.sample}_*_pseudobulk_phenotype.tsv"), emit: "cellranger_pseudobulk_phenotype"

    script:
    """
    pseudobulk.py ${count_table} ${sample_meta.sample} ${airr_schema} --sobject_gex ${sobject_gex}
    """
}