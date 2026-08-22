
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ANNOTATE_PROCESS; ANNOTATE_SORT_CDR3; ANNOTATE_DEDUPLICATE_CDR3_TRBV; ANNOTATE_DEDUPLICATE_CDR3 } from '../../modules/local/annotate'
include { OLGA_CONCATENATE as ANNOTATE_OLGA_CONCATENATE; OLGA_CALCULATE as ANNOTATE_OLGA_CALCULATE} from '../../modules/local/olga'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANNOTATE_INGEST  (bulk-only entrypoint)

    Runs ANNOTATE_PROCESS per raw input sample and builds the concatenated CDR3 table +
    per-sample pre-filter-stats sidecar that ANNOTATE (below) needs. Called only from
    TCRTOOLKIT_BULK, on samples fresh out of INPUT_CHECK/CONVERT.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ANNOTATE_INGEST {
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

    concat_cdr3 = processed_samples
        .map { _meta, f -> f }
        .collectFile(name: 'concat_cdr3.tsv', keepHeader: true, skip: 1)

    emit:
    processed_samples
    per_sample_stats
    concat_cdr3
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ANNOTATE  (shared entrypoint — both bulk and single-cell feed into this)

    Runs the sort -> dedup -> OLGA chain on an already-concatenated CDR3 table. Takes a
    sample_map already in the canonical clonotype schema (either ANNOTATE_INGEST's output,
    for bulk, or a pseudobulk bridge's output, for single-cell) plus its matching
    pre-concatenated CDR3 table, and emits the same channels regardless of caller. Column
    schema is verified identical across all producers (ANNOTATE_PROCESS, sc_to_cdr3.py,
    vdj_to_bulk.py): junction_aa, v_call, d_call, j_call, duplicate_count,
    junction_aa_length, duplicate_frequency_percent, sequence, sequence_id, junction, sample.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ANNOTATE {
    take:
    sample_map    // channel: [meta, file] already in canonical CDR3 schema
    concat_cdr3   // path:    pre-concatenated CDR3 table

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

    emit:
    processed_samples = sample_map
    concat_cdr3_sorted
    cdr3_pgen = ANNOTATE_OLGA_CONCATENATE.out.cdr3_pgen
    olga_stats = ANNOTATE_OLGA_CONCATENATE.out.cdr3_pgen_stats
        .map { f ->
            def _m = f.readLines()
                .collect{ stats -> stats.split('\t') }
                .collectEntries{ stats -> [(stats[0]): stats[1]] }
        }
        .first()
}
