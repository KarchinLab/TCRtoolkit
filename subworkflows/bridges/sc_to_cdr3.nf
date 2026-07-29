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

    // Build sample_map from unit_map.csv so meta carries patient (+ phenotype).
    // patient_id lets the shared PATIENT step pool per-patient for clustering.
    SC_TO_CDR3.out.unit_map
        .splitCsv(header: true)
        .map { row ->
            def meta = [
                sample    : row.sample,
                patient_id: (row.patient ?: row.sample),
                phenotype : (row.phenotype ?: '')
            ]
            [ meta, file(row.file) ]
        }
        .set { sample_map }

    emit:
    sample_map
    concat_cdr3 = SC_TO_CDR3.out.concat_cdr3
}
