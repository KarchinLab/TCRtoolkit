
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
    processed_samples
        .map { meta, file -> [ meta.patient, file ] }
        .groupTuple()
        .set { patient_groups }

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

    if(params.use_gliph2) {
        GLIPH2_TURBOGLIPH(
            PATIENT_CONCATENATE.out.patient_cdr3
        )
    }
}