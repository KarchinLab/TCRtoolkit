#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { VDJ_TO_BULK } from '../../modules/bridges/vdj_to_bulk.nf'

/*
 * VDJ_TO_BULK_SW
 *
 * Takes VDJ_QC's contigs_after_qc.tsv and produces a TCRtoolkit-compatible
 * sample_map channel (one [meta, file] tuple per sample) plus a synthetic
 * samplesheet file.
 *
 * Used in VDJ-only mode when no GEX Seurat object is provided and
 * TCELL_INTEGRATION is skipped.
 */
workflow VDJ_TO_BULK_SW {
    take:
    contigs_after_qc    // path: contigs_after_qc.tsv from VDJ_QC

    main:
    def ss = params.sample_sheet ? file(params.sample_sheet) : file("${projectDir}/assets/NO_FILE")
    VDJ_TO_BULK(
        contigs_after_qc,
        params.vdj_meta_sample_col ?: 'sample',
        ss
    )

    samplesheet_utf8 = VDJ_TO_BULK.out.samplesheet

    // Parse synthetic samplesheet → Nextflow sample_map channel.
    // row.file is the path as seen INSIDE the task, which only exists afterwards when the
    // task ran in the shared work dir. Under process.scratch (AWS Batch / Cirro) it was a
    // throwaway /tmp/nxf.* dir, so take the file from the declared output and use the CSV
    // only for metadata, pairing the two by file name.
    def tsv_by_name = VDJ_TO_BULK.out.bulk_tsv_files
        .flatten()
        .map { f -> [f.name, f] }

    samplesheet_utf8
        .splitCsv(header: true, sep: ',')
        .map { row -> [row.file.tokenize('/').last(), row.findAll { k, _v -> k != 'file' }] }
        .join(tsv_by_name, failOnMismatch: true)
        .map { _name, meta, file_obj -> [meta, file_obj] }
        .set { sample_map }

    emit:
    sample_map
    samplesheet_utf8
}
