
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

include { INPUT_CHECK }         from '../subworkflows/local/input_check'
include { CONVERT }             from '../subworkflows/local/convert'
include { RESOLVE_SAMPLESHEET } from '../subworkflows/local/resolve_samplesheet'
include { SAMPLE }              from '../subworkflows/local/sample'
include { COMPARE }             from '../subworkflows/local/compare'
include { VALIDATE_PARAMS }     from '../subworkflows/local/validate_params'
include { ANNOTATE }            from '../subworkflows/local/annotate'

include { PSEUDOBULK_PHENOTYPE }from '../subworkflows/local/pseudobulk_phenotype'

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

    if (input_format == 'adaptive') {
        CONVERT(
            INPUT_CHECK.out.sample_map,
            input_format
        )
        .sample_map_converted
        .set { sample_map_final }

    } else if (input_format == 'cellranger') {
        CONVERT(
            INPUT_CHECK.out.sample_map,
            input_format
        )
        .sample_map_converted
        .set { sample_map_final }

        if (params.sobject_gex) {
            // Current SCRATCH-annotate gex input:
            // data/SCRATCH_ANNOTATION:SCTYPE_STATE_ANNOTATION/data/project_T_Cells_annotation_object.RDS
            PSEUDOBULK_PHENOTYPE(
                CONVERT.out.pseudobulk_phenotype_files,
                INPUT_CHECK.out.samplesheet_utf8,
                levels
            )
        }

    } else {
        INPUT_CHECK.out.sample_map
            .set { sample_map_final }
    }

    // --- Main Analysis ---

    RESOLVE_SAMPLESHEET( ch_samplesheet_utf8,
        sample_map_final )

    ANNOTATE( RESOLVE_SAMPLESHEET.out.samplesheet_resolved,
        RESOLVE_SAMPLESHEET.out.all_sample_files )

    // Running sample level analysis
    if (levels.contains('sample')) {
        SAMPLE( sample_map_final ),
        ANNOTATE.out.cdr3_pgen
    }

    // Running comparison analysis
    if (levels.contains('compare')) {
        COMPARE( RESOLVE_SAMPLESHEET.out.samplesheet_resolved,
            RESOLVE_SAMPLESHEET.out.all_sample_files,
            ANNOTATE.out.concat_cdr3_sorted,
            ANNOTATE.out.cdr3_pgen)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
