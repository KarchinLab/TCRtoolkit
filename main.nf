#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { BULK_WORKFLOW }        from './workflows/bulk.nf'
include { SINGLECELL_WORKFLOW }  from './workflows/singlecell.nf'
include { COMBINED_WORKFLOW }    from './workflows/combined.nf'

def resolveMode() {
    if (params.mode) return params.mode.toLowerCase()
    def hasBulk = params.samplesheet != null
    // VDJ contigs alone are sufficient for singlecell mode (VDJ-only path)
    def hasSC   = params.input_vdj_contigs != null
    def hasGEX  = params.input_annotated_object != null
    if (hasBulk && (hasSC || hasGEX)) return 'combined'
    if (hasSC || hasGEX)              return 'singlecell'
    if (hasBulk)                      return 'bulk'
    error """
    No valid inputs detected. Please provide:
      --samplesheet                        (bulk mode)
      --input_vdj_contigs + --sample_sheet (singlecell mode; add --input_annotated_object for full SC mode)
      both of the above                    (combined mode)
    Or set --mode explicitly: bulk | singlecell | combined
    """
}

workflow {
    def mode = resolveMode()

    log.info """
    ╔══════════════════════════════════════════════════════════════╗
    ║                         TCR-Toolkit                          ║
    ║         Single-Cell TCR Repertoire Analysis Pipeline         ║
    ╚══════════════════════════════════════════════════════════════╝
    Mode           : ${mode}
    Project        : ${params.project_name}
    Output dir     : ${params.outdir}
    """

    switch (mode) {
        case 'bulk':
            BULK_WORKFLOW()
            break
        case 'singlecell':
            SINGLECELL_WORKFLOW()
            break
        case 'combined':
            COMBINED_WORKFLOW()
            break
        default:
            error "Unknown --mode '${mode}'. Valid options: bulk | singlecell | combined"
    }
}

workflow.onComplete {
    def status = workflow.success ? 'SUCCESS' : 'FAILED'
    log.info """
    ══════════════════════════════════════════════
    Pipeline ${status}
    Results : ${params.outdir}
    Duration: ${workflow.duration}
    ══════════════════════════════════════════════
    """
}
