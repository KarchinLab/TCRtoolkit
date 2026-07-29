
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CONVERT_ADAPTIVE }                from '../../modules/local/convert/convert_adaptive'
// NOTE: cellranger pseudobulk (PSEUDOBULK_CELLRANGER / PSEUDOBULK_PHENOTYPE_CELLRANGER) was
// relocated to the single-cell modality (SC_TO_CDR3 / VDJ_TO_BULK). Bulk keeps airr + adaptive.

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
    pseudobulk_phenotype_files = channel.empty() // retained (empty) for output compatibility

    if (input_format == 'adaptive') {
        CONVERT_ADAPTIVE(
            sample_map,
            file(params.airr_schema),
            file(params.imgt_lookup)
        )
        sample_map_converted = CONVERT_ADAPTIVE.out.adaptive_convert

    } else {
        sample_map_converted = sample_map
    }

    emit:
    sample_map_converted
    pseudobulk_phenotype_files
}