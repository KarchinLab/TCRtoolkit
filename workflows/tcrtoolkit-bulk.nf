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

include { MAP_PHENOTYPES }      from '../subworkflows/local/map_phenotypes'

include { RESOLVE_SAMPLESHEET as RESOLVE_SAMPLESHEET_PHENO } from '../subworkflows/local/resolve_samplesheet'
include { SAMPLE as SAMPLE_PHENO } from '../subworkflows/local/sample'
include { COMPARE as COMPARE_PHENO } from '../subworkflows/local/compare'

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

    // ---  ---

    // Empty channels used for phenotype analysis
    ch_phenotype_files_transformed = Channel.empty() // Initialize empty transformed channel
    ch_phenotype_samplesheet = Channel.empty() // Initialize empty samplesheet for pheno files

    if (input_format in ['adaptive', 'cellranger']) {
        AIRR_CONVERT(
            INPUT_CHECK.out.sample_map,
            input_format,
            params.sobject_gex
        )
        
        sample_map_final = AIRR_CONVERT.out.sample_map_converted

        // --- Phenotype file handling ---
        MAP_PHENOTYPES(
            AIRR_CONVERT.out.pseudobulk_phenotype_files,
            INPUT_CHECK.out.samplesheet_utf8
        )

        ch_phenotype_files_transformed = MAP_PHENOTYPES.out.files_transformed
        ch_phenotype_samplesheet = MAP_PHENOTYPES.out.samplesheet_pheno

    } else {
        sample_map_final = INPUT_CHECK.out.sample_map
    }

    // --- Main Analysis ---
    RESOLVE_SAMPLESHEET( 
        INPUT_CHECK.out.samplesheet_utf8,
        sample_map_final 
    )

    if (levels.contains('sample')) {
        SAMPLE( sample_map_final )
    }

    if (levels.contains('compare')) {
        COMPARE( 
            RESOLVE_SAMPLESHEET.out.samplesheet_resolved,
            RESOLVE_SAMPLESHEET.out.all_sample_files
        )
    }
    
    // --- Phenotype Analysis ---
    
    // These processes will be skipped if their input channels are empty
    RESOLVE_SAMPLESHEET_PHENO(
        ch_phenotype_samplesheet,
        ch_phenotype_files_transformed
    )

    if (levels.contains('sample')) {
        SAMPLE_PHENO( ch_phenotype_files_transformed )
    }

    if (levels.contains('compare')) {
        COMPARE_PHENO(
            RESOLVE_SAMPLESHEET_PHENO.out.samplesheet_resolved,
            RESOLVE_SAMPLESHEET_PHENO.out.all_sample_files
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