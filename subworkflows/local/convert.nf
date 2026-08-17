
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CONVERT_ADAPTIVE } from '../../modules/local/convert/convert_adaptive'

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
}
