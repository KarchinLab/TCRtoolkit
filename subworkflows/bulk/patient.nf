include { PATIENT_CONCATENATE; PATIENT_CALC } from '../../modules/bulk/patient/main'
include { GLIPH2_TURBOGLIPH }               from '../../modules/bulk/compare/gliph2'
include { GIANA_CALC }                      from '../../modules/bulk/compare/giana'

workflow PATIENT {
    take:
    processed_samples

    main:
    // Support both bulk samplesheets (patient column) and SC-derived samplesheets
    // (params.patient_col, e.g. 'patient_id').  Fall back to 'patient' if unset.
    def patient_key = params.patient_col ?: 'patient'

    processed_samples
        .map { meta, file -> [ meta[patient_key] ?: meta.patient, file ] }
        .groupTuple()
        .set { patient_groups }

    PATIENT_CONCATENATE( patient_groups )
    PATIENT_CALC( PATIENT_CONCATENATE.out.patient_cdr3 )

    GIANA_CALC(
        PATIENT_CONCATENATE.out.patient_cdr3,
        params.threshold,
        params.threshold_score,
        params.threshold_vgene
    )

    gliph2_details = Channel.empty()
    if (params.use_gliph2) {
        GLIPH2_TURBOGLIPH( PATIENT_CONCATENATE.out.patient_cdr3 )
        gliph2_details = GLIPH2_TURBOGLIPH.out.cluster_member_details
    }

    emit:
    // GIANA_CALC has two positional outputs: [0] VgeneScores, [1] giana.txt
    giana_clusters       = GIANA_CALC.out[1]
    gliph2_cluster_details = gliph2_details
}
