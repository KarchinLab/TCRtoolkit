
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

include { MAP_PHENOTYPES }      from '../subworkflows/local/map_phenotypes'
include { RESOLVE_SAMPLESHEET_PHENO } from '../subworkflows/local/resolve_samplesheet_pheno'
include { SAMPLE_PHENO } from '../subworkflows/local/sample_pheno'
include { COMPARE_PHENO } from '../subworkflows/local/compare_pheno'

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
    ch_samplesheet_utf8 = INPUT_CHECK.out.samplesheet_utf8

    // Initialize empty channels
    ch_phenotype_files_transformed = Channel.empty() 
    ch_phenotype_samplesheet = Channel.empty() 

    if (input_format in ['adaptive', 'cellranger']) {
        def sobject_gex_file = params.sobject_gex ? file(params.sobject_gex) : []

        AIRR_CONVERT( INPUT_CHECK.out.sample_map,
            input_format, sobject_gex_file
            )
            .sample_map_converted
            .set { sample_map_final }

        if (params.sobject_gex) {
            MAP_PHENOTYPES(
                AIRR_CONVERT.out.pseudobulk_phenotype_files,
                ch_samplesheet_utf8
            )
            ch_phenotype_files_transformed = MAP_PHENOTYPES.out.files_transformed
            ch_phenotype_samplesheet = MAP_PHENOTYPES.out.samplesheet_pheno
        }

    } else {
        INPUT_CHECK.out.sample_map
            .set { sample_map_final }
    }

    // --- Main Analysis ---

    RESOLVE_SAMPLESHEET( ch_samplesheet_utf8,
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

    // --- Phenotype Analysis ---

    if (params.sobject_gex) {

        RESOLVE_SAMPLESHEET_PHENO(
            ch_phenotype_files_transformed
        )

        if (levels.contains('sample')) {
            SAMPLE_PHENO( ch_phenotype_files_transformed, ch_phenotype_samplesheet ) // ch_phenotype_samplesheet // RESOLVE_SAMPLESHEET_PHENO.out.samplesheet_resolved
        }

        if (levels.contains('compare')) {
            COMPARE_PHENO(
                ch_phenotype_samplesheet,
                RESOLVE_SAMPLESHEET_PHENO.out.all_sample_files
            )
        }
    }
    
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
