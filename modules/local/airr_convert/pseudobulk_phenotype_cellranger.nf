process PSEUDOBULK_CELLRANGER_PHENOTYPE {
    tag "${sample_meta.sample}"
    label 'process_low'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    input:
    tuple val(sample_meta), path(count_table)
    path airr_schema
    path sobject_gex

    output:
    tuple val(sample_meta), path("${sample_meta.sample}_*_pseudobulk_phenotype.tsv"), emit: "cellranger_pseudobulk_phenotype"

    script:
    // Dynamically add phenotype arguments only if sobject_gex is provided (i.e., not null)
    def phenotype_args = sobject_gex ? "--phenotype --sobject_gex ${sobject_gex}" : ""

    """
    pseudobulk.py \\
        ${count_table} \\
        ${sample_meta.sample} \\
        ${airr_schema} \\
        ${phenotype_args}
    """
}