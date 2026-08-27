
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INPUT_CHECK }         from '../subworkflows/local/input_check'
include { CONVERT }             from '../subworkflows/local/convert'
include { VALIDATE_PARAMS }     from '../subworkflows/local/validate_params'
include { ANNOTATE_INGEST }     from '../subworkflows/local/annotate'
include { BULKTCR_ANALYSIS }    from '../subworkflows/local/bulktcr_analysis'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow TCRTOOLKIT_BULK {
    VALIDATE_PARAMS()

    println("Running TCRTOOLKIT_BULK workflow...")

    // Construct levels list from the run_sample, run_compare, and run_patient parameters
    def levels = []
    if (params.run_sample) levels << 'sample'
    if (params.run_compare) levels << 'compare'
    if (params.run_patient) levels << 'patient'
    if (params.run_convert) levels << 'convert'

    def input_format = params.input_format.toLowerCase()

    // Validate
    if (params.workflow_level) {
        println("[WARN] workflow_level is set for deprecation, please use run_<level> flags. Overriding flags with workflow_level.")
        levels = params.workflow_level.toLowerCase().tokenize(',')
    }

    if (levels.contains('convert') && input_format != 'adaptive') {
        println("\u001B[33m[WARN]\u001B[0m To run Convert workflow, please specify a valid convertible --input_format (adaptive)")
        if (!levels.contains('sample') && !levels.contains('compare')) {
            return
        }
    }

    if (levels.contains('patient')) {
        def samplesheet_header = file(params.samplesheet).readLines().first().split(',')
        def has_patient = samplesheet_header.contains('patient')

        if (!has_patient) {
            println("\u001B[33m[WARN]\u001B[0m Patient workflow was specified but metadata was not found in samplesheet; please specify patient IDs for samples using the 'patient' column or remove 'patient' from workflow_level.")
            return
        }
    }

    // Checking input tables
    INPUT_CHECK( file(params.samplesheet) )

    if (input_format == 'adaptive') {
        CONVERT(INPUT_CHECK.out.sample_map, input_format)
        sample_map_final = CONVERT.out
    } else {
        sample_map_final = INPUT_CHECK.out.sample_map
    }

    // template_discovery_brief.qmd stages AIRR-converted files only when CONVERT ran
    // (adaptive); template_discovery_brief.qmd's VDJdb section otherwise reads the raw
    // input directly, which already has AIRR-standard frequency columns.
    def convert_files = (input_format == 'adaptive')
        ? CONVERT.out.map { _meta, f -> f }.collect()
        : channel.value([])

    // Bulk reports are sample-centric. Compare- and patient-dependent sections are
    // added only when those workflow levels are present. Change this once reports do
    // not always require sample workflow.
    def run_reports = levels.contains('sample')

    // --- Main Analysis ---
    if (levels.intersect(['sample','patient','compare'])) {
        ANNOTATE_INGEST( sample_map_final )

        BULKTCR_ANALYSIS(
            ANNOTATE_INGEST.out.processed_samples,
            ANNOTATE_INGEST.out.per_sample_stats,
            ANNOTATE_INGEST.out.concat_cdr3,
            levels,
            run_reports,
            convert_files
        )
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
