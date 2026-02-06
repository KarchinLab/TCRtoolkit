
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ANNOTATE_CONCATENATE; ANNOTATE_SORT_CDR3; ANNOTATE_DEDUPLICATE_CDR3_TRBV; ANNOTATE_DEDUPLICATE_CDR3 } from '../../modules/local/annotate'
include { OLGA_CONCATENATE; OLGA_CALCULATE } from '../../modules/local/olga'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ANNOTATE {
    take:
    samplesheet_resolved
    all_sample_files

    main:
    ANNOTATE_CONCATENATE( samplesheet_resolved,
        all_sample_files )

    ANNOTATE_SORT_CDR3( ANNOTATE_CONCATENATE.out.concat_cdr3 )
    concat_cdr3_sorted = ANNOTATE_SORT_CDR3.out.concat_cdr3_sorted

    ANNOTATE_DEDUPLICATE_CDR3_TRBV( concat_cdr3_sorted )

    ANNOTATE_DEDUPLICATE_CDR3(
        ANNOTATE_DEDUPLICATE_CDR3_TRBV.out.unique_cdr3_trbv
    )

    OLGA_CALCULATE(
        ANNOTATE_DEDUPLICATE_CDR3.out.unique_cdr3
            .splitText(by: params.olga_chunk_length, file: true)
    )

    OLGA_CONCATENATE (
        OLGA_CALCULATE.out.pgen_chunk
            .collectFile(
                name: 'olga_pgen_body.tsv',
                sort: { f ->
                    def m = (f.name =~ /\.(\d+)\.txt$/)
                    m ? m[0][1].toInteger() : 0
                }
            )
    )
    cdr3_pgen = OLGA_CONCATENATE.out.cdr3_pgen

    emit:
    concat_cdr3_sorted
    cdr3_pgen
}