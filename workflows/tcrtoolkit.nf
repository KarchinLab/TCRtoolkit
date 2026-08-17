
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


workflow TCRTOOLKIT {
    VALIDATE_PARAMS()

    println("Running TCRTOOLKIT workflow...")

    // Split the workflow_level parameter into a list of levels
    def levels = params.workflow_level.toLowerCase().tokenize(',')
    def input_format = params.input_format.toLowerCase()

    // Validate
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
        sample_map_final = CONVERT.out.sample_map_converted
    } else {
        sample_map_final = INPUT_CHECK.out.sample_map
    }

    // template_discovery_brief.qmd stages AIRR-converted files only when CONVERT ran
    // (adaptive); template_discovery_brief.qmd's VDJdb section otherwise reads the raw
    // input directly, which already has AIRR-standard frequency columns.
    def convert_files = (input_format == 'adaptive')
        ? CONVERT.out.sample_map_converted.map { _meta, f -> f }.collect()
        : channel.value([])

    // --- Main Analysis ---
    if (levels.intersect(['sample','patient','compare'])) {
        ANNOTATE_INGEST( sample_map_final )

        BULKTCR_ANALYSIS(
            ANNOTATE_INGEST.out.processed_samples,
            ANNOTATE_INGEST.out.per_sample_stats,
            ANNOTATE_INGEST.out.concat_cdr3,
            levels,
            true,
            convert_files
        )
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
