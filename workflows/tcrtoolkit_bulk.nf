
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
include { CLONOTYPE_CLUSTERING } from '../subworkflows/local/clonotype_clustering'
include { BULK_TO_EXPORT }      from '../modules/bridges/bulk_to_export'
include { MASTER_SUMMARY_SW }   from '../subworkflows/scratch/master_summary'
include { REPORT }              from '../subworkflows/local/report'

def enabled(x) { x == null || x == true }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow TCRTOOLKIT_BULK {
    VALIDATE_PARAMS()

    println("Running TCRTOOLKIT_BULK workflow...")

    // .toBoolean(): CLI overrides can arrive as the string "false", which is truthy in Groovy.
    def levels = []
    if (params.run_sample?.toString()?.toBoolean() ?: false) levels << 'sample'
    if (params.run_compare?.toString()?.toBoolean() ?: false) levels << 'compare'
    if (params.run_patient?.toString()?.toBoolean() ?: false) levels << 'patient'
    if (params.run_convert?.toString()?.toBoolean() ?: false) levels << 'convert'

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

    // [meta, file] pairs for template_discovery_brief.qmd's VDJdb section.
    def convert_files = (input_format == 'adaptive')
        ? CONVERT.out.collect(flat: false)
        : INPUT_CHECK.out.sample_map.collect(flat: false)

    // Bulk reports are sample-centric. Compare- and patient-dependent sections are
    // added only when those workflow levels are present. Change this once reports do
    // not always require sample workflow.
    def run_reports = levels.contains('sample') && enabled(params.run_reports)

    // --- Main Analysis ---
    if (levels.intersect(['sample','patient','compare'])) {
        ANNOTATE_INGEST( sample_map_final )

        BULKTCR_ANALYSIS(
            ANNOTATE_INGEST.out.processed_samples,
            ANNOTATE_INGEST.out.per_sample_stats,
            ANNOTATE_INGEST.out.concat_cdr3,
            levels,
            run_reports,
            convert_files,
            channel.value(file(params.samplesheet)),
            channel.value(file("${projectDir}/assets/NO_FILE")),
            true
        )

        // Report notebooks; NO_FILE = no single-cell data, so the SC sections stay empty.
        REPORT(
            BULKTCR_ANALYSIS.out.reports,
            channel.value(file(params.samplesheet)),
            channel.value(file("${projectDir}/assets/NO_FILE"))
        )

        // --- Clonotype export -> cluster mapping -> consensus -> Master Summary ---
        // The same downstream chain the single-cell modality runs, after the cohort
        // comparison. There is no GEX object, so it runs on a clonotype-level export
        // (one row per sample x clonotype; bulk counts are reads/templates, too many to
        // expand into one row each).
        // Needs SAMPLE (TCRdist3) and/or PATIENT (GIANA, GLIPH2) for anything to cluster.
        def has_clusters = levels.contains('sample') || levels.contains('patient')
        if (levels.contains('compare') && has_clusters && enabled(params.run_clonotype_clustering)) {
            def nofile = file("${projectDir}/assets/NO_FILE")
            def ch_project_name = channel.value(params.project_name)

            BULK_TO_EXPORT(
                BULKTCR_ANALYSIS.out.concat_cdr3_sorted,
                file(params.samplesheet),
                'clonotype'
            )

            CLONOTYPE_CLUSTERING(
                channel.value(nofile),
                false,
                BULK_TO_EXPORT.out.export_cells,
                BULKTCR_ANALYSIS.out.giana_clusters,
                BULKTCR_ANALYSIS.out.gliph2_cluster_details,
                BULKTCR_ANALYSIS.out.tcrdist_clone_df,
                BULKTCR_ANALYSIS.out.tcrdist_output,
                ch_project_name
            )

            if (enabled(params.run_master_summary)) {
                MASTER_SUMMARY_SW(
                    channel.value(nofile),
                    CLONOTYPE_CLUSTERING.out.rep_export,
                    channel.value([]),
                    channel.value([]),
                    channel.empty()
                        .mix(BULKTCR_ANALYSIS.out.sample_stats, BULKTCR_ANALYSIS.out.v_family, BULKTCR_ANALYSIS.out.j_family)
                        .collect().ifEmpty([]),
                    BULKTCR_ANALYSIS.out.shared_cdr3.flatten().collect().ifEmpty([]),
                    CLONOTYPE_CLUSTERING.out.rollup_tables,
                    BULKTCR_ANALYSIS.out.repertoire_tables.flatten().collect().ifEmpty([]),
                    channel.value([]),
                    channel.value([]),
                    CLONOTYPE_CLUSTERING.out.consensus_tables.flatten().collect().ifEmpty([]),
                    channel.value([]),
                    channel.value([]),
                    CLONOTYPE_CLUSTERING.out.consensus_report.collect().ifEmpty([nofile]),
                    ch_project_name
                )
            }
        }
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
