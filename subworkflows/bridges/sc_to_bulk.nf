#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SC_TO_BULK } from '../../modules/bridges/sc_to_bulk.nf'

/*
 * SC_TO_BULK_SW
 *
 * Takes the per-cell export_cells.tsv from SCRATCH T-cell integration and
 * produces a TCRtoolkit-compatible sample_map channel (one [meta, file] tuple
 * per sample) plus a synthetic samplesheet file for report modules.
 *
 * The metadata columns carried over (patient, condition, timepoint, batch)
 * are read from the first cell per sample in export_cells.tsv.
 */
workflow SC_TO_BULK_SW {
    take:
    export_cells    // path: per-cell TSV from TCELL_INTEGRATION

    main:
    def meta_cols = [
        params.patient_col,
        params.condition_col,
        params.timepoint_col,
        params.batch_col
    ].findAll { it }.join(',')

    // TCELL_INTEGRATION always exports a standardized 'sample' column;
    // params.sample_col is the Seurat metadata name, not the export column name.
    SC_TO_BULK(
        export_cells,
        'sample',
        meta_cols
    )

    samplesheet_utf8 = SC_TO_BULK.out.samplesheet

    // Parse synthetic samplesheet → Nextflow sample_map channel
    samplesheet_utf8
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta     = row.findAll { k, _v -> k != 'file' }
            def file_obj = file(row.file)
            return [meta, file_obj]
        }
        .set { sample_map }

    emit:
    sample_map
    samplesheet_utf8
}
