#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SC_TO_CDR3 } from '../../modules/bridges/sc_to_cdr3.nf'

/*
 * SC_TO_CDR3_SW
 *
 * Wraps SC_TO_CDR3 and emits:
 *   sample_map   — channel of [meta, file] tuples (one per sample)
 *   concat_cdr3  — single path to concatenated CDR3 file
 *
 * Used in full-SC mode as a replacement for SC_TO_BULK_SW + ANNOTATE_PROCESS.
 */
workflow SC_TO_CDR3_SW {
    take:
    export_cells    // path: tcr_export_cells_with_embedding.tsv from TCELL_INTEGRATION

    main:
    SC_TO_CDR3( export_cells )

    // Build sample_map channel: [meta, file] per sample
    SC_TO_CDR3.out.per_sample
        .flatten()
        .map { f ->
            def sample = f.name.replaceAll('_cdr3\\.tsv$', '')
            [[sample: sample], f]
        }
        .set { sample_map }

    emit:
    sample_map
    concat_cdr3 = SC_TO_CDR3.out.concat_cdr3
}
