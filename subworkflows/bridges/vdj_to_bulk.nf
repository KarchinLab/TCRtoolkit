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
