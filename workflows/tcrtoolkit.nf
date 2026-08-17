
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
include { ANNOTATE_INGEST; ANNOTATE } from '../subworkflows/local/annotate'
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
        ANNOTATE_INGEST( sample_map_final )
        ANNOTATE( ANNOTATE_INGEST.out.processed_samples, ANNOTATE_INGEST.out.concat_cdr3 )
    }

    // Running sample level analysis
    if (levels.contains('sample')) {
        SAMPLE(
            ANNOTATE.out.processed_samples,
            ANNOTATE_INGEST.out.per_sample_stats,
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

    // Report - works on channel of tuples [notebook template, files to stage, staged directory layout]
    ch_reports = channel.empty()

    def sample_stats_agg = SAMPLE.out.sample_csv
        .collectFile(name: "sample_stats.csv", keepHeader: true, skip: 1, sort: true)

    // Only stage AIRR-converted files when CONVERT ran (adaptive/cellranger);
    // template_discovery_brief.qmd's VDJdb section otherwise reads the raw
    // input directly, which already has AIRR-standard frequency columns.
    def convert_files = (input_format == 'adaptive' || input_format == 'cellranger')
        ? CONVERT.out.sample_map_converted.map { _meta, f -> f }.collect()
        : channel.value([])

    // template_discovery_brief.qmd includes a single generic template_pheno.qmd -
    // RENDER_NOTEBOOK stages whichever real notebook applies under that shared
    // name (see staged_layout below), since Quarto's {{< include >}} shortcode is
    // a static textual splice with no native runtime if/else. pheno_sc needs the
    // per-cell/per-phenotype pseudobulk files, which only exist for cellranger
    // input with sobject_gex supplied; everything else gets pheno_bulk, which
    // only needs TCRPHENO output (produced for every sample regardless of format).
    def use_pheno_sc = (input_format == 'cellranger' && params.sobject_gex)
    def pheno_notebook = use_pheno_sc ? file(params.template_pheno_sc) : file(params.template_pheno_bulk)
    // .collect() on a channel that never emits (channel.empty(), which is what
    // CONVERT.out.pseudobulk_phenotype_files is outside the cellranger+sobject_gex
    // case) never emits either - not even an empty list - which would silently
    // starve every downstream .combine() in ch_discovery_report and make
    // RENDER_NOTEBOOK(template_discovery_brief) never get scheduled at all (no
    // error, just silently skipped - confirmed directly with a minimal repro).
    // Guard with the same channel.value([]) fallback already used for convert_files.
    def pseudobulk_pheno_files = use_pheno_sc
        ? CONVERT.out.pseudobulk_phenotype_files.map { _meta, files -> files }.flatten().collect()
        : channel.value([])

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

    // Discovery brief reads its inputs from a project_dir/<project_name>/
    // <subdir>/<file> layout, so staged_layout tells RENDER_NOTEBOOK where to
    // symlink each staged file - each entry is a [dest_path, source_basename]
    // pair (source and dest basenames usually match, but not always - see
    // gliph2 below).
    //
    // .collect()-produced list channels are wrapped via `.map { l -> [l] }`
    // before every .combine() below - otherwise .combine() flattens the
    // list's contents into the tuple instead of keeping it as one element.
    //
    // .combine() (no `by:`) is a safe pairing here, not a risky cross-product:
    // every channel below is a whole-run aggregate that emits exactly one item
    // per pipeline run (a .collectFile() result or a .collect()'d list), so
    // there's nothing to key by - there's only ever one item on each side.
    ch_discovery_report = sample_stats_agg
        .combine(ANNOTATE.out.concat_cdr3_sorted)
        .combine(COMPARE.out.shared_cdr3)
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
            // template_details_part1.qmd pulls in template_sample.qmd via a
            // {{< include >}} shortcode resolved relative to its own directory,
            // so that sibling file has to be staged alongside it too - Nextflow
            // only stages the single notebook file named in the tuple otherwise.
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

    // template_giana.qmd/template_gliph.qmd need patient-level clustering
    // outputs that only exist when 'patient' is in workflow_level (GIANA's
    // concat crashes outright otherwise), so template_patient_clustering.qmd
    // is resolved to an "on" (includes both) or "off" (placeholder) wrapper,
    // mirroring the template_pheno.qmd trick above.
    def details_part2_base = sample_stats_agg
        .combine(ANNOTATE.out.concat_cdr3_sorted)
        .combine(COMPARE.out.shared_cdr3)

    def run_patient_clustering = levels.contains('patient')
    def patient_clustering_notebook = run_patient_clustering
        ? file(params.template_patient_clustering_on)
        : file(params.template_patient_clustering_off)

    // template_details_part2.qmd pulls in these files via {{< include >}}
    // shortcodes resolved relative to its own directory, so they have to be
    // staged alongside it too - see the equivalent part1 comment above.
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
                // gliph2 files keep their patient-prefixed basename (e.g.
                // "patientA_all_motifs.txt") all the way through, since
                // template_gliph.qmd reads that same name from each patient's
                // subdir - source and dest basenames match here.
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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
