/*
 * BULK workflow (Scenario 1)
 * Runs the full TCRtoolkit bulk repertoire analysis.
 * Input: AIRR / Adaptive / CellRanger TSV files via --samplesheet
 */

include { INPUT_CHECK }          from '../subworkflows/bulk/input_check.nf'
include { CONVERT }              from '../subworkflows/bulk/convert.nf'
include { ANNOTATE }             from '../subworkflows/bulk/annotate.nf'
include { SAMPLE }               from '../subworkflows/bulk/sample.nf'
include { PATIENT }              from '../subworkflows/bulk/patient.nf'
include { COMPARE }              from '../subworkflows/bulk/compare.nf'
include { PSEUDOBULK_PHENOTYPE } from '../subworkflows/bulk/pseudobulk_phenotype.nf'
include { VALIDATE_PARAMS }      from '../subworkflows/bulk/validate_params.nf'

workflow BULK_WORKFLOW {

    VALIDATE_PARAMS()

    def levels       = params.workflow_level.toLowerCase().tokenize(',')
    def input_format = params.input_format.toLowerCase()

    if (levels.contains('convert') && !['adaptive', 'cellranger'].contains(input_format)) {
        log.warn "[WARN] To run Convert workflow, specify a valid --input_format (adaptive or cellranger)"
        if (!levels.contains('sample') && !levels.contains('compare')) return
    }

    def runPatient = levels.contains('patient')
    if (runPatient) {
        def header = file(params.samplesheet).readLines().first().split(',')
        if (!header.contains('patient')) {
            log.warn "[WARN] 'patient' level requested but no 'patient' column found in samplesheet — skipping patient level."
            runPatient = false
        }
    }

    INPUT_CHECK( file(params.samplesheet) )

    if (input_format == 'adaptive' || input_format == 'cellranger') {
        CONVERT(
            INPUT_CHECK.out.sample_map,
            input_format
        )
        sample_map_final = CONVERT.out.sample_map_converted

        if (input_format == 'cellranger' && params.sobject_gex) {
            PSEUDOBULK_PHENOTYPE(
                CONVERT.out.pseudobulk_phenotype_files,
                INPUT_CHECK.out.samplesheet_utf8,
                levels
            )
        }
    } else {
        sample_map_final = INPUT_CHECK.out.sample_map
    }

    if (levels.intersect(['sample', 'patient', 'compare'])) {
        ANNOTATE( sample_map_final )
    }

    if (levels.contains('sample')) {
        SAMPLE(
            sample_map_final,
            ANNOTATE.out.cdr3_pgen,
            ANNOTATE.out.olga_stats
        )
    }

    if (runPatient) {
        PATIENT( ANNOTATE.out.processed_samples )
    }

    if (levels.contains('compare')) {
        COMPARE(
            ANNOTATE.out.concat_cdr3_sorted,
            ANNOTATE.out.cdr3_pgen
        )
    }
}
