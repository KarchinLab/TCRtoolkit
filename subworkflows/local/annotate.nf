
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ANNOTATE_PROCESS; ANNOTATE_SORT_CDR3; ANNOTATE_DEDUPLICATE_CDR3_TRBV; ANNOTATE_DEDUPLICATE_CDR3 } from '../../modules/local/annotate'
include { OLGA_CONCATENATE as ANNOTATE_OLGA_CONCATENATE; OLGA_CALCULATE as ANNOTATE_OLGA_CALCULATE} from '../../modules/local/olga'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ANNOTATE {
    take:
    sample_map

    main:
    ANNOTATE_PROCESS( sample_map )
    
    processed_samples = ANNOTATE_PROCESS.out.process
    per_sample_stats = ANNOTATE_PROCESS.out.pre_filter_stats

    per_sample_stats
        .map { _meta, stats_file -> stats_file }
        .collectFile(
            name: "pre_filter_stats.csv",
            keepHeader: true,
            skip: 1,
            sort: true,
            storeDir: "${params.outdir}/sample"
        )
    
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

    ANNOTATE_OLGA_CONCATENATE (
        ANNOTATE_OLGA_CALCULATE.out.pgen_chunk
            .collectFile(
                name: 'olga_pgen_body.tsv',
                sort: { f ->
                    def m = (f.name =~ /\.(\d+)\.txt$/)
                    m ? m[0][1].toInteger() : 0
                }
            )
    )
    cdr3_pgen = ANNOTATE_OLGA_CONCATENATE.out.cdr3_pgen
    olga_stats = ANNOTATE_OLGA_CONCATENATE.out.cdr3_pgen_stats
        .map { f ->
            def _m = f.readLines()
                .collect{ stats -> stats.split('\t') }
                .collectEntries{ stats -> [(stats[0]): stats[1]] }
        }
        .first()

    emit:
    processed_samples
    per_sample_stats
    concat_cdr3_sorted
    cdr3_pgen
    olga_stats
}