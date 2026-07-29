
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PATIENT_CONCATENATE; PATIENT_CALC } from '../../modules/local/patient'
include { GLIPH2_TURBOGLIPH } from '../../modules/local/compare/gliph2'
include { GIANA_CALC    } from '../../modules/local/compare/giana'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PATIENT {
    take:
    processed_samples

    main:
    // Grouping key is configurable for the single-cell modality (e.g. 'patient_id').
    // When params.patient_col is unset (all bulk runs), this resolves to meta.patient —
    // identical to the previous behavior.
    def patient_key = params.patient_col ?: 'patient'
    def patient_groups = processed_samples
        .map { meta, file -> [ meta[patient_key] ?: meta.patient, file ] }
        .groupTuple()

    PATIENT_CONCATENATE ( patient_groups )

    PATIENT_CALC( PATIENT_CONCATENATE.out.patient_cdr3 )

    // TODO: disabling plotting until notebook updated
    // COMPARE_PLOT( samplesheet_resolved,
    //             COMPARE_CALC.out.jaccard_mat,
    //             COMPARE_CALC.out.sorensen_mat,
    //             COMPARE_CALC.out.morisita_mat,
    //             file(params.compare_stats_template),
    //             params.project_name,
    //             all_sample_files
    //             )

    GIANA_CALC(
        PATIENT_CONCATENATE.out.patient_cdr3,
        params.threshold,
        params.threshold_score,
        params.threshold_vgene
    )

    // Capture GLIPH2 cluster details into a defaulted channel so the emit below is always
    // valid; when use_gliph2 is false it stays empty. Bulk mode ignores these emits.
    def gliph2_details = Channel.empty()
    if(params.use_gliph2) {
        GLIPH2_TURBOGLIPH(
            PATIENT_CONCATENATE.out.patient_cdr3
        )
        gliph2_details = GLIPH2_TURBOGLIPH.out.cluster_member_details_named
    }

    emit:
    // Additive outputs consumed only by the single-cell modality (CLUSTER_TO_SC).
    // GIANA_CALC's second positional output is the giana.txt cluster file.
    giana_clusters         = GIANA_CALC.out[1]
    gliph2_cluster_details = gliph2_details
}