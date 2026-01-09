
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CONVERT_ADAPTIVE } from '../../modules/local/airr_convert/convert_adaptive'
include { PSEUDOBULK_CELLRANGER } from '../../modules/local/airr_convert/pseudobulk_cellranger'
include { PSEUDOBULK_CELLRANGER_PHENOTYPE } from '../../modules/local/airr_convert/pseudobulk_phenotype_cellranger'

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
            .adaptive_convert
            .set { sample_map_converted }
    } else if (input_format == 'cellranger') {
        PSEUDOBULK_CELLRANGER(
            sample_map,
            params.airr_schema
        )
            .cellranger_pseudobulk
            .set { sample_map_converted }
        
        // Pseudobulk by phenotype
        if(sobject_gex) {
            PSEUDOBULK_CELLRANGER_PHENOTYPE(
                sample_map,
                params.airr_schema,
                file(sobject_gex) // Can remove file()
            )
            // Capture phenotype outputs
            ch_pseudobulk_phenotype = PSEUDOBULK_CELLRANGER_PHENOTYPE.out.cellranger_pseudobulk_phenotype
        }
    }

    emit:
    sample_map_converted
    pseudobulk_phenotype_files = ch_pseudobulk_phenotype
}