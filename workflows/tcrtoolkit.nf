
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
include { SAMPLE }              from '../subworkflows/local/sample'
include { PATIENT }             from '../subworkflows/local/patient'
include { COMPARE }             from '../subworkflows/local/compare'
include { VALIDATE_PARAMS }     from '../subworkflows/local/validate_params'
include { ANNOTATE }            from '../subworkflows/local/annotate'
include { REPORT }              from '../subworkflows/local/report'

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

    } else if (input_format == 'cellranger') {
        CONVERT(INPUT_CHECK.out.sample_map, input_format)
        sample_map_final = CONVERT.out.sample_map_converted

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
        sample_map_final = INPUT_CHECK.out.sample_map
    }

    // --- Main Analysis ---
    if (levels.intersect(['sample','patient','compare'])) {
        ANNOTATE( sample_map_final )
    }

    // Running sample level analysis
    if (levels.contains('sample')) {
        SAMPLE(
            ANNOTATE.out.processed_samples,
            ANNOTATE.out.per_sample_stats,
            ANNOTATE.out.cdr3_pgen,
            ANNOTATE.out.olga_stats
        )
    }

    // Running patient analysis
    if (levels.contains('patient')) {
        PATIENT( ANNOTATE.out.processed_samples )
    }

    // Running comparison analysis
    if (levels.contains('compare')) {
        COMPARE(
            ANNOTATE.out.concat_cdr3_sorted,
            ANNOTATE.out.cdr3_pgen
        )
    }

    // Report - works on channel of tuples [report name, [report files]]
    ch_reports = channel.empty()

    // QC report requires sample-level aggregate outputs.
    ch_qc_report = SAMPLE.out.sample_csv
        .collectFile(name: "sample_stats.csv", keepHeader: true, skip: 1, sort: true)
        .combine(ANNOTATE.out.concat_cdr3_sorted)
        .map { sample_stats_csv, concat_cdr3_sorted ->
            tuple(
                file(params.template_qc),
                [sample_stats_csv,
                concat_cdr3_sorted]
            )
        }
    ch_reports = ch_reports.mix(ch_qc_report)

    // Another report
    // ch_reports = ch_reports.mix( ... )

    REPORT( ch_reports )


}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
