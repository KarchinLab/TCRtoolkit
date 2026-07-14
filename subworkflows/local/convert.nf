
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CONVERT_ADAPTIVE }                from '../../modules/local/convert/convert_adaptive'
include { PSEUDOBULK_CELLRANGER }           from '../../modules/local/convert/pseudobulk_cellranger'
include { PSEUDOBULK_PHENOTYPE_CELLRANGER } from '../../modules/local/convert/pseudobulk_phenotype_cellranger'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CONVERT {
    take:
    sample_map
    input_format

    main:    
    pseudobulk_phenotype_files = channel.empty() // Initialize empty channel for phenotype files

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