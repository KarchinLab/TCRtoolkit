
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

    The single shared bulk-TCR analysis engine. Both TCRTOOLKIT (bulk mode) and
    SINGLECELL_WORKFLOW (single-cell mode, on its pseudobulked data) invoke this
    subworkflow rather than importing ANNOTATE/SAMPLE/PATIENT/COMPARE individually.

    `levels` preserves TCRTOOLKIT's --workflow_level partial-run feature (which of
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

    if (levels.contains('sample')) {
        SAMPLE(
            ANNOTATE.out.processed_samples,
            per_sample_stats,
            ANNOTATE.out.cdr3_pgen,
            ANNOTATE.out.olga_stats
        )
    }

    if (levels.contains('patient')) {
        PATIENT( ANNOTATE.out.processed_samples )
    }

    if (levels.contains('compare')) {
        COMPARE(
            ANNOTATE.out.concat_cdr3_sorted,
            ANNOTATE.out.cdr3_pgen
        )
    }

    if (run_reports) {
        ch_reports = channel.empty()

        def sample_stats_agg = SAMPLE.out.sample_csv
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

        ch_discovery_report = sample_stats_agg
            .combine(ANNOTATE.out.concat_cdr3_sorted)
            .combine(COMPARE.out)
            .combine(SAMPLE.out.tcrdist_files.map { l -> [l] })
            .combine(SAMPLE.out.vdjdb_files.map { l -> [l] })
            .combine(convert_files.map { l -> [l] })
            .combine(SAMPLE.out.tcrpheno_files.map { l -> [l] })
            .combine(pseudobulk_pheno_files.map { l -> [l] })
            .map { sample_stats_csv, concat_cdr3_sorted, shared_cdr3, tcrdist_files, vdjdb_files, convert_files_l, tcrpheno_files, pseudobulk_files ->
                def report_files = [sample_stats_csv, concat_cdr3_sorted, shared_cdr3, pheno_notebook] +
                    tcrdist_files + vdjdb_files + convert_files_l + tcrpheno_files + pseudobulk_files
                def staged_layout = [
                    ["${params.project_name}/sample/${sample_stats_csv.name}", sample_stats_csv.name],
                    ["${params.project_name}/annotate/${concat_cdr3_sorted.name}", concat_cdr3_sorted.name],
                    ["${params.project_name}/tcrsharing/${shared_cdr3.name}", shared_cdr3.name],
                    ["template_pheno.qmd", pheno_notebook.name]
                ] + tcrdist_files.collect { f -> ["${params.project_name}/tcrdist3/${f.name}", f.name] } +
                    vdjdb_files.collect { f -> ["${params.project_name}/vdjdb/${f.name}", f.name] } +
                    convert_files_l.collect { f -> ["${params.project_name}/convert/${f.name}", f.name] } +
                    tcrpheno_files.collect { f -> ["${params.project_name}/tcrpheno/${f.name}", f.name] } +
                    pseudobulk_files.collect { f -> ["${params.project_name}/pseudobulk/${f.name}", f.name] }
                tuple(
                    file(params.template_discovery_brief),
                    report_files,
                    staged_layout
                )
            }
        ch_reports = ch_reports.mix(ch_discovery_report)

        ch_details_part1_report = sample_stats_agg
            .combine(ANNOTATE.out.concat_cdr3_sorted)
            .combine(SAMPLE.out.v_family)
            .combine(SAMPLE.out.j_family)
            .combine(SAMPLE.out.tcrdist_files.map { l -> [l] })
            .combine(SAMPLE.out.olga_files.map { l -> [l] })
            .combine(SAMPLE.out.vdjdb_files.map { l -> [l] })
            .combine(SAMPLE.out.convergence_files.map { l -> [l] })
            .map { sample_stats_csv, concat_cdr3_sorted, v_family, j_family, tcrdist_files, olga_files, vdjdb_files, convergence_files ->
                def include_files = [file("${file(params.template_details_part1).parent}/template_sample.qmd")]
                def report_files = [sample_stats_csv, concat_cdr3_sorted, v_family, j_family] +
                    tcrdist_files + olga_files + vdjdb_files + convergence_files + include_files
                def staged_layout = [
                    ["${params.project_name}/sample/${sample_stats_csv.name}", sample_stats_csv.name],
                    ["${params.project_name}/annotate/${concat_cdr3_sorted.name}", concat_cdr3_sorted.name],
                    ["${params.project_name}/sample/${v_family.name}", v_family.name],
                    ["${params.project_name}/sample/${j_family.name}", j_family.name]
                ] + tcrdist_files.collect { f -> ["${params.project_name}/tcrdist3/${f.name}", f.name] } +
                    olga_files.collect { f -> ["${params.project_name}/olga/${f.name}", f.name] } +
                    vdjdb_files.collect { f -> ["${params.project_name}/vdjdb/${f.name}", f.name] } +
                    convergence_files.collect { f -> ["${params.project_name}/convergence/${f.name}", f.name] }
                tuple(
                    file(params.template_details_part1),
                    report_files,
                    staged_layout
                )
            }
        ch_reports = ch_reports.mix(ch_details_part1_report)

        def details_part2_base = sample_stats_agg
            .combine(ANNOTATE.out.concat_cdr3_sorted)
            .combine(COMPARE.out)

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
                .combine(PATIENT.out.giana_files.map { l -> [l] })
                .combine(PATIENT.out.gliph2_all_motifs.map { l -> [l] })
                .combine(PATIENT.out.gliph2_clone_network.map { l -> [l] })
                .combine(PATIENT.out.gliph2_cluster_member_details.map { l -> [l] })
                .combine(PATIENT.out.gliph2_global_similarities.map { l -> [l] })
                .map { sample_stats_csv, concat_cdr3_sorted, shared_cdr3, giana_files,
                       all_motifs_pairs, clone_network_pairs, cluster_member_pairs, global_sim_pairs ->
                    def report_files = [sample_stats_csv, concat_cdr3_sorted, shared_cdr3] + giana_files +
                        all_motifs_pairs.collect { p -> p[1] } +
                        clone_network_pairs.collect { p -> p[1] } +
                        cluster_member_pairs.collect { p -> p[1] } +
                        global_sim_pairs.collect { p -> p[1] } +
                        part2_include_files
                    def staged_layout = [
                        ["${params.project_name}/sample/${sample_stats_csv.name}", sample_stats_csv.name],
                        ["${params.project_name}/annotate/${concat_cdr3_sorted.name}", concat_cdr3_sorted.name],
                        ["${params.project_name}/tcrsharing/${shared_cdr3.name}", shared_cdr3.name]
                    ] + giana_files.collect { f -> ["${params.project_name}/giana/${f.name}", f.name] } +
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
                .map { sample_stats_csv, concat_cdr3_sorted, shared_cdr3 ->
                    tuple(
                        file(params.template_details_part2),
                        [sample_stats_csv, concat_cdr3_sorted, shared_cdr3] + part2_include_files,
                        [
                            ["${params.project_name}/sample/${sample_stats_csv.name}", sample_stats_csv.name],
                            ["${params.project_name}/annotate/${concat_cdr3_sorted.name}", concat_cdr3_sorted.name],
                            ["${params.project_name}/tcrsharing/${shared_cdr3.name}", shared_cdr3.name]
                        ] + part2_staged_layout_extra
                    )
                }
        }
        ch_reports = ch_reports.mix(ch_details_part2_report)

        REPORT( ch_reports )
    }

    emit:
    processed_samples  = ANNOTATE.out.processed_samples
    concat_cdr3_sorted = ANNOTATE.out.concat_cdr3_sorted
    cdr3_pgen          = ANNOTATE.out.cdr3_pgen
    olga_stats         = ANNOTATE.out.olga_stats

    sample_csv        = levels.contains('sample') ? SAMPLE.out.sample_csv        : channel.empty()
    tcrdist_clone_df  = levels.contains('sample') ? SAMPLE.out.tcrdist_clone_df  : channel.empty()
    tcrdist_output    = levels.contains('sample') ? SAMPLE.out.tcrdist_output    : channel.empty()
    v_family          = levels.contains('sample') ? SAMPLE.out.v_family          : channel.empty()
    j_family          = levels.contains('sample') ? SAMPLE.out.j_family          : channel.empty()
    tcrdist_files     = levels.contains('sample') ? SAMPLE.out.tcrdist_files     : channel.value([])
    olga_files        = levels.contains('sample') ? SAMPLE.out.olga_files        : channel.value([])
    vdjdb_files       = levels.contains('sample') ? SAMPLE.out.vdjdb_files       : channel.value([])
    convergence_files = levels.contains('sample') ? SAMPLE.out.convergence_files : channel.value([])
    tcrpheno_files    = levels.contains('sample') ? SAMPLE.out.tcrpheno_files    : channel.value([])

    giana_clusters                = levels.contains('patient') ? PATIENT.out.giana_clusters                : channel.empty()
    gliph2_cluster_details         = levels.contains('patient') ? PATIENT.out.gliph2_cluster_details         : channel.empty()
    giana_files                    = levels.contains('patient') ? PATIENT.out.giana_files                    : channel.value([])
    gliph2_all_motifs              = levels.contains('patient') ? PATIENT.out.gliph2_all_motifs              : channel.value([])
    gliph2_clone_network           = levels.contains('patient') ? PATIENT.out.gliph2_clone_network           : channel.value([])
    gliph2_cluster_member_details  = levels.contains('patient') ? PATIENT.out.gliph2_cluster_member_details  : channel.value([])
    gliph2_global_similarities     = levels.contains('patient') ? PATIENT.out.gliph2_global_similarities     : channel.value([])

    shared_cdr3 = levels.contains('compare') ? COMPARE.out : channel.empty()
}
