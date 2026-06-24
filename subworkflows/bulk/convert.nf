include { CONVERT_ADAPTIVE }                from '../../modules/bulk/convert/convert_adaptive'
include { PSEUDOBULK_CELLRANGER }           from '../../modules/bulk/convert/pseudobulk_cellranger'
include { PSEUDOBULK_PHENOTYPE_CELLRANGER } from '../../modules/bulk/convert/pseudobulk_phenotype_cellranger'

workflow CONVERT {
    take:
    sample_map
    input_format

    main:
    pseudobulk_phenotype_files = Channel.empty()

    if (input_format == 'adaptive') {
        CONVERT_ADAPTIVE(
            sample_map,
            file(params.airr_schema),
            file(params.imgt_lookup)
        )
        sample_map_converted = CONVERT_ADAPTIVE.out.adaptive_convert

    } else if (input_format == 'cellranger') {
        PSEUDOBULK_CELLRANGER(
            sample_map,
            file(params.airr_schema)
        )
        sample_map_converted = PSEUDOBULK_CELLRANGER.out.cellranger_pseudobulk

        if (params.sobject_gex) {
            PSEUDOBULK_PHENOTYPE_CELLRANGER(
                sample_map,
                file(params.airr_schema),
                file(params.sobject_gex)
            )
            pseudobulk_phenotype_files = PSEUDOBULK_PHENOTYPE_CELLRANGER.out.cellranger_pseudobulk_phenotype
        }
    } else {
        sample_map_converted = sample_map
    }

    emit:
    sample_map_converted
    pseudobulk_phenotype_files
}
