
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

    main:    
    if (input_format == 'adaptive') {
        CONVERT_ADAPTIVE(
            sample_map,
            params.airr_schema,
            params.imgt_lookup
        )
            .adaptive_convert
            .set { sample_map_converted }
    } else if (input_format == 'cellranger') {
        PSEUDOBULK_CELLRANGER(
            sample_map,
            params.airr_schema
        )
            .cellranger_pseudobulk
            // .cellranger_pseudobulk_phenotype // Changed
            .set { sample_map_converted }
    }

    emit:
    sample_map_converted
}