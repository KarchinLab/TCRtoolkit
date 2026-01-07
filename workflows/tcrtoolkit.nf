
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { INPUT_CHECK }         from '../subworkflows/local/input_check'
include { AIRR_CONVERT }        from '../subworkflows/local/airr_convert'
include { RESOLVE_SAMPLESHEET } from '../subworkflows/local/resolve_samplesheet'
include { SAMPLE }              from '../subworkflows/local/sample'
include { COMPARE }             from '../subworkflows/local/compare'
include { VALIDATE_PARAMS }     from '../subworkflows/local/validate_params'

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
    if (levels.contains('convert') && !['adaptive', 'cellranger'].contains(input_format)) {
        println("\u001B[33m[WARN]\u001B[0m To run Convert workflow, please specify a valid convertible --input_format (adaptive or cellranger)")
        if (!levels.contains('sample') && !levels.contains('compare')) {
            return
        }
    }

    // Checking input tables
    INPUT_CHECK( file(params.samplesheet) )

    if (input_format in ['adaptive', 'cellranger']) {
        AIRR_CONVERT( INPUT_CHECK.out.sample_map,
            input_format
            )
            .sample_map_converted
            .set { sample_map_final }
    } else {
        INPUT_CHECK.out.sample_map
            .set { sample_map_final }
    }

    RESOLVE_SAMPLESHEET( INPUT_CHECK.out.samplesheet_utf8,
        sample_map_final )

    // Running sample level analysis
    if (levels.contains('sample')) {
        SAMPLE( sample_map_final )
    }

    // Running comparison analysis
    if (levels.contains('compare')) {
        COMPARE( RESOLVE_SAMPLESHEET.out.samplesheet_resolved,
            RESOLVE_SAMPLESHEET.out.all_sample_files)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
