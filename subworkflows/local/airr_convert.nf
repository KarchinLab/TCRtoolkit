
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CONVERT_ADAPTIVE } from '../../modules/local/airr_convert/convert_adaptive'
include { PSEUDOBULK_CELLRANGER } from '../../modules/local/airr_convert/pseudobulk_cellranger'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow AIRR_CONVERT {
    take:
    sample_map
    input_format
    sobject_gex 

    main:
    ch_pseudobulk_phenotype = Channel.empty() // Initialize empty channel for phenotype files

    if (input_format == 'adaptive') {
        CONVERT_ADAPTIVE(
            sample_map,
            params.airr_schema,
            params.imgt_lookup
        )
        sample_map_converted = CONVERT_ADAPTIVE.out.adaptive_convert
    } else if (input_format == 'cellranger') {
        PSEUDOBULK_CELLRANGER(
            sample_map,
            params.airr_schema,
            sobject_gex 
        )
        sample_map_converted = PSEUDOBULK_CELLRANGER.out.cellranger_pseudobulk
        ch_pseudobulk_phenotype = PSEUDOBULK_CELLRANGER.out.cellranger_pseudobulk_phenotype
    }

    emit:
    sample_map_converted
    pseudobulk_phenotype_files = ch_pseudobulk_phenotype
}


// ======================================

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Validate pipeline parameters
def checkPathParamList = [ params.samplesheet, params.sobject_gex ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.samplesheet) { samplesheet = file(params.samplesheet) } else { exit 1, 'Samplesheet not specified. Please, provide a --samplesheet=/path/to/samplesheet.csv !' }
if (params.outdir) { outdir = params.outdir } else { exit 1, 'Output directory not specified. Please, provide a --outdir=/path/to/outdir !' }

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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow TCRTOOLKIT_BULK {

    println("Running TCRTOOLKIT_BULK workflow...")

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

    // MODIFIED: This entire 'if/else' block is updated
    ch_phenotype_files = Channel.empty() // Initialize empty channel

    if (input_format in ['adaptive', 'cellranger']) {
        AIRR_CONVERT(
            INPUT_CHECK.out.sample_map,
            input_format,
            params.sobject_gex // <-- 1. Pass sobject_gex parameter
        )
        
        // 2. Capture outputs using the .out syntax
        sample_map_final = AIRR_CONVERT.out.sample_map_converted
        ch_phenotype_files = AIRR_CONVERT.out.pseudobulk_phenotype_files // Capture the new channel

    } else {
        sample_map_final = INPUT_CHECK.out.sample_map
    }

    RESOLVE_SAMPLESHEET(
        INPUT_CHECK.out.samplesheet_utf8,
        sample_map_final
    )

    // Running sample level analysis
    if (levels.contains('sample')) {
        SAMPLE( sample_map_final )
    }

    // Running comparison analysis
    if (levels.contains('compare')) {
        COMPARE(
            RESOLVE_SAMPLESHEET.out.samplesheet_resolved,
            RESOLVE_SAMPLESHEET.out.all_sample_files
        )
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// workflow.onComplete {

//     log.info(workflow.success ? "Finished tcrtoolkit-bulk!" : "Please check your inputs.")

// }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/