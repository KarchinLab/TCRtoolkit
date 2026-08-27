
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ANNOTATE } from './annotate'
include { SAMPLE }   from './sample'
include { PATIENT }  from './patient'
include { COMPARE }  from './compare'
include { REPORT }   from './report'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    BULKTCR_ANALYSIS

    The single shared bulk-TCR analysis engine. Both TCRTOOLKIT_BULK (bulk mode) and
    TCRTOOLKIT_SC (single-cell mode, on its pseudobulked data) invoke this
    subworkflow rather than importing ANNOTATE/SAMPLE/PATIENT/COMPARE individually.

    `levels` preserves TCRTOOLKIT_BULK's --workflow_level partial-run feature (which of
    sample/patient/compare to run). `run_reports` gates the 4 bulk report notebooks:
    bulk's report templates assume samplesheet columns (origin/timepoint) that
    single-cell-derived samplesheets don't reliably carry, so SC callers pass
    run_reports: false and use their own REPERTOIRE/MASTER_SUMMARY reports instead.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BULKTCR_ANALYSIS {
    take:
    sample_map          // channel [meta, file]: per-unit table, canonical clonotype schema
    per_sample_stats    // channel [meta, file]: pre_filter_stats.csv per unit
    concat_cdr3         // path: pre-concatenated CDR3 table across units
    levels              // list, e.g. ['sample','patient','compare'] - which stages to run
    run_reports         // bool: render the 4 bulk report notebooks?
    convert_files       // channel: list of converted-input files for the discovery-brief
                         //   report (bulk-only; SC callers pass channel.value([]))

    main:
    ANNOTATE( sample_map, concat_cdr3 )

    def ch_sample_csv = channel.empty()
    def ch_tcrdist_clone_df = channel.empty()
    def ch_tcrdist_output = channel.empty()
    def ch_v_family = channel.empty()
    def ch_j_family = channel.empty()
    def ch_tcrdist_files = channel.value([])
    def ch_olga_files = channel.value([])
    def ch_vdjdb_files = channel.value([])
    def ch_convergence_files = channel.value([])
    def ch_tcrpheno_files = channel.value([])

    def ch_giana_clusters = channel.empty()
    def ch_gliph2_cluster_details = channel.empty()
    def ch_giana_files = channel.value([])
    def ch_gliph2_all_motifs = channel.value([])
    def ch_gliph2_clone_network = channel.value([])
    def ch_gliph2_cluster_member_details = channel.value([])
    def ch_gliph2_global_similarities = channel.value([])

    def ch_shared_cdr3 = channel.empty()

    if (levels.contains('sample')) {
        SAMPLE(
            ANNOTATE.out.processed_samples,
            per_sample_stats,
            ANNOTATE.out.cdr3_pgen,
            ANNOTATE.out.olga_stats
        )

        ch_sample_csv = SAMPLE.out.sample_csv
        ch_tcrdist_clone_df = SAMPLE.out.tcrdist_clone_df
        ch_tcrdist_output = SAMPLE.out.tcrdist_output
        ch_v_family = SAMPLE.out.v_family
        ch_j_family = SAMPLE.out.j_family
        ch_tcrdist_files = SAMPLE.out.tcrdist_files
        ch_olga_files = SAMPLE.out.olga_files
        ch_vdjdb_files = SAMPLE.out.vdjdb_files
        ch_convergence_files = SAMPLE.out.convergence_files
        ch_tcrpheno_files = SAMPLE.out.tcrpheno_files
    }

    if (levels.contains('patient')) {
        PATIENT( ANNOTATE.out.processed_samples )

        ch_giana_clusters = PATIENT.out.giana_clusters
        ch_gliph2_cluster_details = PATIENT.out.gliph2_cluster_details
        ch_giana_files = PATIENT.out.giana_files
        ch_gliph2_all_motifs = PATIENT.out.gliph2_all_motifs
        ch_gliph2_clone_network = PATIENT.out.gliph2_clone_network
        ch_gliph2_cluster_member_details = PATIENT.out.gliph2_cluster_member_details
        ch_gliph2_global_similarities = PATIENT.out.gliph2_global_similarities
    }

    if (levels.contains('compare')) {
        COMPARE(
            ANNOTATE.out.concat_cdr3_sorted,
            ANNOTATE.out.cdr3_pgen
        )

        ch_shared_cdr3 = COMPARE.out
    }

    if (run_reports && levels.contains('sample')) {
        ch_reports = channel.empty()

        def sample_stats_agg = ch_sample_csv
            .collectFile(name: "sample_stats.csv", keepHeader: true, skip: 1, sort: true)

        // input_format can now only be adaptive/airr (cellranger removed), so the
        // single-cell phenotype notebook branch is dead: always use the bulk notebook,
        // and there are never pseudobulk-phenotype files to stage.
        def pheno_notebook = file(params.template_pheno_bulk)
        def pseudobulk_pheno_files = channel.value([])

        ch_qc_report = sample_stats_agg
            .combine(ANNOTATE.out.concat_cdr3_sorted)
            .map { sample_stats_csv, concat_cdr3_sorted ->
                tuple(
                    file(params.template_qc),
                    [sample_stats_csv,
                    concat_cdr3_sorted],
                    []
                )
            }
        ch_reports = ch_reports.mix(ch_qc_report)

        ch_details_part1_report = sample_stats_agg
            .combine(ANNOTATE.out.concat_cdr3_sorted)
            .combine(ch_v_family)
            .combine(ch_j_family)
            .combine(ch_tcrdist_files.map { l -> [l] })
            .combine(ch_olga_files.map { l -> [l] })
            .combine(ch_vdjdb_files.map { l -> [l] })
            .combine(ch_convergence_files.map { l -> [l] })
            .map { sample_stats_csv, concat_cdr3_sorted, v_family_file, j_family_file, tcrdist_files_l, olga_files_l, vdjdb_files_l, convergence_files_l ->
                def include_files = [file("${file(params.template_details_part1).parent}/template_sample.qmd")]
                def report_files = [sample_stats_csv, concat_cdr3_sorted, v_family_file, j_family_file] +
                    tcrdist_files_l + olga_files_l + vdjdb_files_l + convergence_files_l + include_files
                def staged_layout = [
                    ["${params.project_name}/sample/${sample_stats_csv.name}", sample_stats_csv.name],
                    ["${params.project_name}/annotate/${concat_cdr3_sorted.name}", concat_cdr3_sorted.name],
                    ["${params.project_name}/sample/${v_family_file.name}", v_family_file.name],
                    ["${params.project_name}/sample/${j_family_file.name}", j_family_file.name]
                ] + tcrdist_files_l.collect { f -> ["${params.project_name}/tcrdist3/${f.name}", f.name] } +
                    olga_files_l.collect { f -> ["${params.project_name}/olga/${f.name}", f.name] } +
                    vdjdb_files_l.collect { f -> ["${params.project_name}/vdjdb/${f.name}", f.name] } +
                    convergence_files_l.collect { f -> ["${params.project_name}/convergence/${f.name}", f.name] }
                tuple(
                    file(params.template_details_part1),
                    report_files,
                    staged_layout
                )
            }
        ch_reports = ch_reports.mix(ch_details_part1_report)

        if (levels.contains('compare')) {
            ch_discovery_report = sample_stats_agg
                .combine(ANNOTATE.out.concat_cdr3_sorted)
                .combine(ch_shared_cdr3)
                .combine(ch_tcrdist_files.map { l -> [l] })
                .combine(ch_vdjdb_files.map { l -> [l] })
                .combine(convert_files.map { l -> [l] })
                .combine(ch_tcrpheno_files.map { l -> [l] })
                .combine(pseudobulk_pheno_files.map { l -> [l] })
                .map { sample_stats_csv, concat_cdr3_sorted, shared_cdr3_file, tcrdist_files_l, vdjdb_files_l, convert_files_l, tcrpheno_files_l, pseudobulk_files_l ->
                    def report_files = [sample_stats_csv, concat_cdr3_sorted, shared_cdr3_file, pheno_notebook] +
                        tcrdist_files_l + vdjdb_files_l + convert_files_l + tcrpheno_files_l + pseudobulk_files_l
                    def staged_layout = [
                        ["${params.project_name}/sample/${sample_stats_csv.name}", sample_stats_csv.name],
                        ["${params.project_name}/annotate/${concat_cdr3_sorted.name}", concat_cdr3_sorted.name],
                        ["${params.project_name}/tcrsharing/${shared_cdr3_file.name}", shared_cdr3_file.name],
                        ["template_pheno.qmd", pheno_notebook.name]
                    ] + tcrdist_files_l.collect { f -> ["${params.project_name}/tcrdist3/${f.name}", f.name] } +
                        vdjdb_files_l.collect { f -> ["${params.project_name}/vdjdb/${f.name}", f.name] } +
                        convert_files_l.collect { f -> ["${params.project_name}/convert/${f.name}", f.name] } +
                        tcrpheno_files_l.collect { f -> ["${params.project_name}/tcrpheno/${f.name}", f.name] } +
                        pseudobulk_files_l.collect { f -> ["${params.project_name}/pseudobulk/${f.name}", f.name] }
                    tuple(
                        file(params.template_discovery_brief),
                        report_files,
                        staged_layout
                    )
                }
            ch_reports = ch_reports.mix(ch_discovery_report)

            def details_part2_base = sample_stats_agg
                .combine(ANNOTATE.out.concat_cdr3_sorted)
                .combine(ch_shared_cdr3)

            def run_patient_clustering = levels.contains('patient')
            def patient_clustering_notebook = run_patient_clustering
                ? file(params.template_patient_clustering_on)
                : file(params.template_patient_clustering_off)

            def part2_notebooks_dir = file(params.template_details_part2).parent
            def part2_include_files = [
                file("${part2_notebooks_dir}/template_overlap.qmd"),
                file("${part2_notebooks_dir}/template_sharing.qmd"),
                file("${part2_notebooks_dir}/template_giana.qmd"),
                file("${part2_notebooks_dir}/template_gliph.qmd"),
                patient_clustering_notebook
            ]
            def part2_staged_layout_extra = [
                ["template_patient_clustering.qmd", patient_clustering_notebook.name]
            ]

            if (run_patient_clustering) {
                ch_details_part2_report = details_part2_base
                    .combine(ch_giana_files.map { l -> [l] })
                    .combine(ch_gliph2_all_motifs.map { l -> [l] })
                    .combine(ch_gliph2_clone_network.map { l -> [l] })
                    .combine(ch_gliph2_cluster_member_details.map { l -> [l] })
                    .combine(ch_gliph2_global_similarities.map { l -> [l] })
                    .map { sample_stats_csv, concat_cdr3_sorted, shared_cdr3_file, giana_files_l,
                           all_motifs_pairs, clone_network_pairs, cluster_member_pairs, global_sim_pairs ->
                        def report_files = [sample_stats_csv, concat_cdr3_sorted, shared_cdr3_file] + giana_files_l +
                            all_motifs_pairs.collect { p -> p[1] } +
                            clone_network_pairs.collect { p -> p[1] } +
                            cluster_member_pairs.collect { p -> p[1] } +
                            global_sim_pairs.collect { p -> p[1] } +
                            part2_include_files
                        def staged_layout = [
                            ["${params.project_name}/sample/${sample_stats_csv.name}", sample_stats_csv.name],
                            ["${params.project_name}/annotate/${concat_cdr3_sorted.name}", concat_cdr3_sorted.name],
                            ["${params.project_name}/tcrsharing/${shared_cdr3_file.name}", shared_cdr3_file.name]
                        ] + giana_files_l.collect { f -> ["${params.project_name}/giana/${f.name}", f.name] } +
                            all_motifs_pairs.collect { p -> ["${params.project_name}/gliph2/${p[0]}/${p[1].name}", p[1].name] } +
                            clone_network_pairs.collect { p -> ["${params.project_name}/gliph2/${p[0]}/${p[1].name}", p[1].name] } +
                            cluster_member_pairs.collect { p -> ["${params.project_name}/gliph2/${p[0]}/${p[1].name}", p[1].name] } +
                            global_sim_pairs.collect { p -> ["${params.project_name}/gliph2/${p[0]}/${p[1].name}", p[1].name] } +
                            part2_staged_layout_extra
                        tuple(
                            file(params.template_details_part2),
                            report_files,
                            staged_layout
                        )
                    }
            } else {
                ch_details_part2_report = details_part2_base
                    .map { sample_stats_csv, concat_cdr3_sorted, shared_cdr3_file ->
                        tuple(
                            file(params.template_details_part2),
                            [sample_stats_csv, concat_cdr3_sorted, shared_cdr3_file] + part2_include_files,
                            [
                                ["${params.project_name}/sample/${sample_stats_csv.name}", sample_stats_csv.name],
                                ["${params.project_name}/annotate/${concat_cdr3_sorted.name}", concat_cdr3_sorted.name],
                                ["${params.project_name}/tcrsharing/${shared_cdr3_file.name}", shared_cdr3_file.name]
                            ] + part2_staged_layout_extra
                        )
                    }
            }
            ch_reports = ch_reports.mix(ch_details_part2_report)
        }

        REPORT( ch_reports )
    }

    emit:
    processed_samples  = ANNOTATE.out.processed_samples
    concat_cdr3_sorted = ANNOTATE.out.concat_cdr3_sorted
    cdr3_pgen          = ANNOTATE.out.cdr3_pgen
    olga_stats         = ANNOTATE.out.olga_stats

    sample_csv        = ch_sample_csv
    tcrdist_clone_df  = ch_tcrdist_clone_df
    tcrdist_output    = ch_tcrdist_output
    v_family          = ch_v_family
    j_family          = ch_j_family
    tcrdist_files     = ch_tcrdist_files
    olga_files        = ch_olga_files
    vdjdb_files       = ch_vdjdb_files
    convergence_files = ch_convergence_files
    tcrpheno_files    = ch_tcrpheno_files

    giana_clusters                = ch_giana_clusters
    gliph2_cluster_details        = ch_gliph2_cluster_details
    giana_files                   = ch_giana_files
    gliph2_all_motifs             = ch_gliph2_all_motifs
    gliph2_clone_network          = ch_gliph2_clone_network
    gliph2_cluster_member_details = ch_gliph2_cluster_member_details
    gliph2_global_similarities    = ch_gliph2_global_similarities

    shared_cdr3 = ch_shared_cdr3
}
