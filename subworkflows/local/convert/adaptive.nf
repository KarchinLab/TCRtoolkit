include { CONVERT_ADAPTIVE } from '../../../modules/local/convert/convert_adaptive'

workflow ADAPTIVE {
    take:
    sample_map

    main:    
    CONVERT_ADAPTIVE(
        sample_map,
        params.airr_schema,
        params.imgt_lookup
    )

    emit:
    CONVERT_ADAPTIVE.out.adaptive_convert
}