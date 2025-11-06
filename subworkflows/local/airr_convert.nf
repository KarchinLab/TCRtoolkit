/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CONVERT_ADAPTIVE } from '../../modules/local/airr_convert/convert_adaptive'
include { PSEUDOBULK_CELLRANGER } from '../../modules/local/airr_convert/pseudobulk_cellranger'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow AIRR_CONVERT {
    take:
    sample_map
    input_format
    sobject_gex 

    main:
    ch_pseudobulk_phenotype = Channel.empty() // Initialize empty channel for phenotype files

    if (input_format == 'adaptive') {
        CONVERT_ADAPTIVE(
            sample_map,
            params.airr_schema,
            params.imgt_lookup
        )
        sample_map_converted = CONVERT_ADAPTIVE.out.adaptive_convert
    } else if (input_format == 'cellranger') {
        PSEUDOBULK_CELLRANGER(
            sample_map,
            params.airr_schema,
            sobject_gex
        )
        // Capture outputs 
        sample_map_converted = PSEUDOBULK_CELLRANGER.out.cellranger_pseudobulk
        ch_pseudobulk_phenotype = PSEUDOBULK_CELLRANGER.out.cellranger_pseudobulk_phenotype
    }

    emit:
    sample_map_converted
    pseudobulk_phenotype_files = ch_pseudobulk_phenotype
}