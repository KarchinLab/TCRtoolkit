
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PATIENT_CONCATENATE } from '../../modules/local/patient'
include { GLIPH2_TURBOGLIPH } from '../../modules/local/compare/gliph2'
include { GIANA_CALC    } from '../../modules/local/compare/giana'
include { TCRDIST3_MATRIX } from '../../modules/local/patient/tcrdist3'

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

    // TCRdist3 on all of a patient's clones at once (the scope GIANA / GLIPH2 use), so
    // similar TCRs are linked across that patient's samples.
    def tcrdist_groups = processed_samples
        .map { meta, file -> [ meta[patient_key] ?: meta.patient, meta.sample, file ] }
        .groupTuple()
    TCRDIST3_MATRIX(
        tcrdist_groups,
        params.matrix_sparsity,
        params.distance_metric,
        file(params.db_path)
    )

    // Within-patient overlap (Jaccard / Sorensen / Morisita) is the patient pairs of
    // REPERTOIRE_METRICS' overlap_samples.tsv.

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
            PATIENT_CONCATENATE.out.patient_cdr3,
            params.local_min_pvalue,
            params.simulation_depth,
            params.kmer_min_depth,
            params.local_min_OVE,
            params.all_aa_interchangeable
        )
    }

    emit:
    // Patient matrices and clone tables, for CLUSTER_TO_SC's TCRdist3 clusters.
    tcrdist_patient_clones = TCRDIST3_MATRIX.out.patient_clones
    tcrdist_patient_matrix = TCRDIST3_MATRIX.out.patient_matrix
    // Per-sample blocks plus the patient files, staged into the report notebooks' tcrdist3/.
    tcrdist_files          = TCRDIST3_MATRIX.out.sample_matrices
                                 .mix( TCRDIST3_MATRIX.out.sample_clone_dfs,
                                       TCRDIST3_MATRIX.out.patient_clones,
                                       TCRDIST3_MATRIX.out.patient_matrix.map { _patient, f -> f } )
                                 .flatten()
                                 .collect()

    // Additive outputs consumed only by the single-cell modality (CLUSTER_TO_SC).
    // GIANA_CALC's second positional output is the giana.txt cluster file.
    giana_clusters         = GIANA_CALC.out[1]
    gliph2_cluster_details = params.use_gliph2 ? GLIPH2_TURBOGLIPH.out.cluster_member_details_named : channel.empty()

    // .collect() flattens tuple(val, path) emissions by default (e.g.
    // [patientA, fileA, patientB, fileB] instead of [[patientA, fileA], ...]),
    // which corrupts the [patient, file] pair indexing used downstream in
    // workflows/tcrtoolkit_bulk.nf - flat: false preserves the pair shape.
    giana_files                   = GIANA_CALC.out.giana_output.collect()
    gliph2_all_motifs             = params.use_gliph2 ? GLIPH2_TURBOGLIPH.out.all_motifs.collect(flat: false) : channel.value([])
    gliph2_clone_network          = params.use_gliph2 ? GLIPH2_TURBOGLIPH.out.clone_network.collect(flat: false) : channel.value([])
    gliph2_cluster_member_details = params.use_gliph2 ? GLIPH2_TURBOGLIPH.out.cluster_member_details.collect(flat: false) : channel.value([])
    gliph2_global_similarities    = params.use_gliph2 ? GLIPH2_TURBOGLIPH.out.global_similarities.collect(flat: false) : channel.value([])
}