include { PSEUDOBULK_CELLRANGER } from '../../../modules/local/convert/pseudobulk_cellranger'

workflow PSEUDOBULK {
    take:
    sample_map

    main:    
    PSEUDOBULK_CELLRANGER(
        sample_map,
        params.airr_schema,
    )

    emit:
    PSEUDOBULK_CELLRANGER.out.cellranger_pseudobulk
}