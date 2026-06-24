include { ANNOTATE_PROCESS; ANNOTATE_SORT_CDR3; ANNOTATE_DEDUPLICATE_CDR3_TRBV; ANNOTATE_DEDUPLICATE_CDR3 } from '../../modules/bulk/annotate/main'
include { OLGA_CONCATENATE as ANNOTATE_OLGA_CONCATENATE; OLGA_CALCULATE as ANNOTATE_OLGA_CALCULATE } from '../../modules/bulk/olga/main'

/*
 * ANNOTATE_FROM_CONCAT
 *
 * Variant of ANNOTATE that skips ANNOTATE_PROCESS.
 * Used in full-SC mode where TCELL_INTEGRATION already provides
 * CDR3b/TRBV/TRBJ/counts/sample data via SC_TO_CDR3_SW.
 *
 * Takes:
 *   sample_map   — channel of [meta, cdr3_tsv] tuples (already in final format)
 *   concat_cdr3  — pre-built concatenated CDR3 file
 */
workflow ANNOTATE_FROM_CONCAT {
    take:
    sample_map   // channel: [meta, file] already in CDR3b/TRBV/TRBJ/counts/sample format
    concat_cdr3  // path: concatenated CDR3 file from SC_TO_CDR3

    main:
    ANNOTATE_SORT_CDR3( concat_cdr3 )
    concat_cdr3_sorted = ANNOTATE_SORT_CDR3.out.concat_cdr3_sorted

    ANNOTATE_DEDUPLICATE_CDR3_TRBV( concat_cdr3_sorted )

    ANNOTATE_DEDUPLICATE_CDR3(
        ANNOTATE_DEDUPLICATE_CDR3_TRBV.out.unique_cdr3_trbv
    )

    ANNOTATE_OLGA_CALCULATE(
        ANNOTATE_DEDUPLICATE_CDR3.out.unique_cdr3
            .splitText(by: params.olga_chunk_length, file: true)
    )

    ANNOTATE_OLGA_CONCATENATE(
        ANNOTATE_OLGA_CALCULATE.out.pgen_chunk
            .collectFile(
                name: 'olga_pgen_body.tsv',
                sort: { f ->
                    def m = (f.name =~ /\.(\d+)\.txt$/)
                    m ? m[0][1].toInteger() : 0
                }
            )
    )

    cdr3_pgen  = ANNOTATE_OLGA_CONCATENATE.out.cdr3_pgen
    olga_stats = ANNOTATE_OLGA_CONCATENATE.out.cdr3_pgen_stats
        .map { f ->
            f.readLines()
                .collect { l -> l.split('\t') }
                .collectEntries { l -> [(l[0]): l[1]] }
        }
        .first()

    emit:
    processed_samples  = sample_map
    concat_cdr3_sorted
    cdr3_pgen
    olga_stats
}

workflow ANNOTATE {
    take:
    sample_map

    main:
    ANNOTATE_PROCESS( sample_map )

    processed_samples = ANNOTATE_PROCESS.out.process

    processed_samples
        .map { _meta, file -> file }
        .collectFile(name: 'concat_cdr3.tsv', keepHeader: true, skip: 1)
        .set { concat_cdr3 }

    ANNOTATE_SORT_CDR3( concat_cdr3 )
    concat_cdr3_sorted = ANNOTATE_SORT_CDR3.out.concat_cdr3_sorted

    ANNOTATE_DEDUPLICATE_CDR3_TRBV( concat_cdr3_sorted )

    ANNOTATE_DEDUPLICATE_CDR3(
        ANNOTATE_DEDUPLICATE_CDR3_TRBV.out.unique_cdr3_trbv
    )

    ANNOTATE_OLGA_CALCULATE(
        ANNOTATE_DEDUPLICATE_CDR3.out.unique_cdr3
            .splitText(by: params.olga_chunk_length, file: true)
    )

    ANNOTATE_OLGA_CONCATENATE(
        ANNOTATE_OLGA_CALCULATE.out.pgen_chunk
            .collectFile(
                name: 'olga_pgen_body.tsv',
                sort: { f ->
                    def m = (f.name =~ /\.(\d+)\.txt$/)
                    m ? m[0][1].toInteger() : 0
                }
            )
    )

    cdr3_pgen  = ANNOTATE_OLGA_CONCATENATE.out.cdr3_pgen
    olga_stats = ANNOTATE_OLGA_CONCATENATE.out.cdr3_pgen_stats
        .map { f ->
            f.readLines()
                .collect { l -> l.split('\t') }
                .collectEntries { l -> [(l[0]): l[1]] }
        }
        .first()

    emit:
    processed_samples
    concat_cdr3_sorted
    cdr3_pgen
    olga_stats
}
