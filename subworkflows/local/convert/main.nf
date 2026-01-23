
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ADAPTIVE } from './adaptive'
include { PSEUDOBULK } from './pseudobulk_cellranger'
include { PSEUDOBULK_PHENOTYPE } from './pseudobulk_phenotype'

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
    ch_pseudobulk_phenotype = channel.empty() // Initialize empty channel for phenotype files

    if (input_format == 'adaptive') {
        ADAPTIVE(
            sample_map,
        )
            .sample_map_converted
            .set { sample_map_converted }
    } else if (input_format == 'cellranger') {
        PSEUDOBULK(
            sample_map,
        )
            .sample_map_converted
            .set { sample_map_converted }
        
        // Pseudobulk by phenotype
        if(params.sobject_gex) {
            PSEUDOBULK_PHENOTYPE(
                sample_map,
                file(params.sobject_gex)
            )
            // Capture phenotype outputs            
            ch_pseudobulk_phenotype = PSEUDOBULK_PHENOTYPE.out.ch_pseudobulk_phenotype
        }
    }

    emit:
    sample_map_converted
    pseudobulk_phenotype_files = ch_pseudobulk_phenotype
}