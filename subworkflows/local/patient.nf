
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

    // Each gliph2_* emit is a list of [patient, file] pairs - kept separate per
    // output type (rather than mixed together) so downstream staging can map
    // each pair back to its known target leaf-name (e.g. "all_motifs.txt")
    // without having to parse it back out of the patient-prefixed filename.
    //
    // NOTE: this Nextflow version's strict-syntax workflow output collection
    // requires emit values to be direct .out expressions - referencing a
    // pre-computed local `def` variable in `emit:` fails at definition time
    // with "Missing workflow output parameter", even outside any conditional.
    // So the params.use_gliph2 ternary has to live in the emit line itself,
    // referencing GLIPH2_TURBOGLIPH.out directly.
    if (params.use_gliph2) {
        GLIPH2_TURBOGLIPH(
            PATIENT_CONCATENATE.out.patient_cdr3
        )
    }

    emit:
    // Additive outputs consumed only by the single-cell modality (CLUSTER_TO_SC).
    // GIANA_CALC's second positional output is the giana.txt cluster file.
    giana_clusters         = GIANA_CALC.out[1]
    gliph2_cluster_details = params.use_gliph2 ? GLIPH2_TURBOGLIPH.out.cluster_member_details_named : channel.empty()

    // .collect() flattens tuple(val, path) emissions by default (e.g.
    // [patientA, fileA, patientB, fileB] instead of [[patientA, fileA], ...]),
    // which corrupts the [patient, file] pair indexing used downstream in
    // workflows/tcrtoolkit.nf - flat: false preserves the pair shape.
    giana_files                   = GIANA_CALC.out.giana_output.collect()
    gliph2_all_motifs             = params.use_gliph2 ? GLIPH2_TURBOGLIPH.out.all_motifs.collect(flat: false) : channel.value([])
    gliph2_clone_network          = params.use_gliph2 ? GLIPH2_TURBOGLIPH.out.clone_network.collect(flat: false) : channel.value([])
    gliph2_cluster_member_details = params.use_gliph2 ? GLIPH2_TURBOGLIPH.out.cluster_member_details.collect(flat: false) : channel.value([])
    gliph2_global_similarities    = params.use_gliph2 ? GLIPH2_TURBOGLIPH.out.global_similarities.collect(flat: false) : channel.value([])
}