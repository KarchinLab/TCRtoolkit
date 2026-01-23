include { PSEUDOBULK_PHENOTYPE_CELLRANGER } from '../../../modules/local/convert/pseudobulk_phenotype_cellranger'

workflow PSEUDOBULK_PHENOTYPE {
    take:
    sample_map
    sobject_gex

    main:    
    PSEUDOBULK_PHENOTYPE_CELLRANGER(
        sample_map,
        params.airr_schema,
        sobject_gex
    )
    ch_pseudobulk_phenotype = PSEUDOBULK_PHENOTYPE_CELLRANGER.out.cellranger_pseudobulk_phenotype

    emit:
    ch_pseudobulk_phenotype
}