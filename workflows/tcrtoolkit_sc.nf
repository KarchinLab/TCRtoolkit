/*
 * TCRTOOLKIT_SC — single-cell TCR modality (integrated).
 *
 * Unified spine (see IMPLEMENTATION_SPEC.md §1.3). Both routes converge after their
 * pseudobulk step into: PSEUDOBULK_QC → ANNOTATE → full shared bulk engine
 * (SAMPLE + REPERTOIRE_METRICS + PATIENT + COMPARE + report notebooks) →
 * CLONOTYPE_CLUSTERING (shared with bulk) → summary. REPERTOIRE_METRICS counts paired
 * clones (paired_clones.tsv from the route's bridge) rather than beta CDR3s.
 *
 *   Full SC   (--input_annotated_object provided):
 *     VDJ_QC → TCELL_INTEGRATION → SC_TO_CDR3 (phenotype pseudobulk)
 *            → PSEUDOBULK_QC → ANNOTATE
 *            → SAMPLE + PATIENT + COMPARE + report notebooks   (shared engine)
 *            → CONGA                                            (GEX-only)
 *            → CLUSTER_TO_SC → CONSENSUS → CLUSTER_ROLLUP       (shared, on the Seurat;
 *                                                                consensus adds CoNGA)
 *            → TCRI                                             (GEX-only)
 *            → MASTER_SUMMARY (full)
 *
 *   VDJ-only  (--input_annotated_object absent):
 *     VDJ_QC → VDJ_TO_BULK (pseudobulk from contigs)
 *            → PSEUDOBULK_QC → ANNOTATE
 *            → SAMPLE + PATIENT + COMPARE + report notebooks   (shared engine)
 *            → BULK_TO_EXPORT (per-cell export from clonotypes)
 *            → CLUSTER_TO_SC → CONSENSUS → CLUSTER_ROLLUP       (shared, table-only)
 *            → MASTER_SUMMARY (CoNGA-excluded)
 *            (only TCELL_INTEGRATION / CONGA / TCRI are skipped — they need the GEX substrate)
 */

// ── Single-cell subworkflows (SCRATCH) ────────────────────────────────────
include { VDJ_QC_SW }            from '../subworkflows/scratch/vdj_qc.nf'
include { TCELL_INTEGRATION_SW } from '../subworkflows/scratch/tcell_integration.nf'
include { CONGA_SW }             from '../subworkflows/scratch/conga.nf'
include { TCRI_SW }             from '../subworkflows/scratch/tcri.nf'
include { MASTER_SUMMARY_SW }    from '../subworkflows/scratch/master_summary.nf'

// ── Bridges (SC ↔ bulk-engine schema conversion) ──────────────────────────
include { SC_TO_CDR3_SW }        from '../subworkflows/bridges/sc_to_cdr3.nf'
include { VDJ_TO_BULK_SW }       from '../subworkflows/bridges/vdj_to_bulk.nf'
include { SC_SAMPLE_STATS }      from '../modules/bridges/sc_sample_stats.nf'
include { SC_REPORT_SAMPLESHEET } from '../modules/bridges/sc_report_samplesheet.nf'
include { REPORT_BUNDLE }        from '../modules/bridges/report_bundle.nf'
include { BULK_TO_EXPORT }       from '../modules/bridges/bulk_to_export.nf'
include { MERGE_VDJ_OBJECT }     from '../modules/bridges/merge_vdj_object.nf'

// ── Shared bulk-TCR analysis engine (main/local — unchanged behavior) ─────
include { PSEUDOBULK_QC_SW }  from '../subworkflows/local/pseudobulk_qc.nf'
include { BULKTCR_ANALYSIS }  from '../subworkflows/local/bulktcr_analysis.nf'
include { CLONOTYPE_CLUSTERING } from '../subworkflows/local/clonotype_clustering.nf'
include { REPORT }            from '../subworkflows/local/report.nf'

// A workflow-local `def enabled = { x -> ... }` closure isn't visible from
// inside nested if-blocks under this Nextflow version's strict-syntax parser
// ("`enabled` is not defined") - a plain top-level function is.
def enabled(x) { x == null || x == true }

// Pull a named file out of a module's collected tables/ channel, falling back to NO_FILE
// so a route where that module did not run still supplies something the report can read
// as "absent" (0 bytes).
def pickFile(ch, fname, nofile) {
    ch.flatten().filter { it.name == fname }.first().ifEmpty(nofile)
}

// Same with a second source. concat() preserves order, so .first() deterministically
// prefers the primary file and falls back otherwise - unlike mix(), which would race.
def pickFileOr(ch, fname, alt, nofile) {
    ch.flatten().filter { it.name == fname }.concat(alt).first().ifEmpty(nofile)
}

workflow TCRTOOLKIT_SC {

    def nofile  = file("${projectDir}/assets/NO_FILE")

    // ── Mandatory inputs (both routes) ────────────────────────────────────
    if (!params.input_vdj_contigs) error "Please provide --input_vdj_contigs"
    if (!params.sample_sheet)      error "Please provide --sample_sheet"

    def vdjOnly = (!params.input_annotated_object ||
                    params.input_annotated_object == '' ||
                    params.input_annotated_object.endsWith('NO_FILE'))

    if (vdjOnly) {
        log.info "==> VDJ-only mode: no GEX annotated object provided. " +
                 "TCELL_INTEGRATION / CONGA / TCRI will be skipped."
    }

    ch_sample_sheet  = channel.fromPath(params.sample_sheet, checkIfExists: true)
    ch_project_name  = channel.value(params.project_name)
    ch_annotated_obj = vdjOnly
        ? channel.fromPath("${projectDir}/assets/NO_FILE")
        : channel.fromPath(params.input_annotated_object, checkIfExists: true)

    // ── Step 1: VDJ QC (both routes) ──────────────────────────────────────
    vdj_qc_out = VDJ_QC_SW( ch_sample_sheet, ch_project_name, ch_annotated_obj )

    // VDJ QC tables/figures kept for the Master Summary (full-SC route)
    def vdj_qc_per_sample_compact         = vdj_qc_out.qc_tables.flatten().filter { f -> f.name == 'vdj_qc_per_sample_compact.tsv'      }.ifEmpty(nofile)
    def vdj_qc_before_after_summary       = vdj_qc_out.qc_tables.flatten().filter { f -> f.name == 'qc_contigs_before_after_summary.tsv' }.ifEmpty(nofile)
    def vdj_qc_sample_sheet_resolved      = vdj_qc_out.qc_tables.flatten().filter { f -> f.name == 'sample_sheet_resolved.tsv'           }.ifEmpty(nofile)
    def vdj_qc_clone_rank_abundance       = vdj_qc_out.qc_tables.flatten().filter { f -> f.name == 'clone_rank_abundance.tsv'            }.ifEmpty(nofile)
    def vdj_qc_before_after_retention_fig = vdj_qc_out.qc_figures.flatten().filter { f -> f.name == 'qc_before_after_retention.png'      }.ifEmpty(nofile)
    def vdj_qc_pairing_bar_fig            = vdj_qc_out.qc_figures.flatten().filter { f -> f.name == 'pairing_bar_by_sample.png'          }.ifEmpty(nofile)
    def vdj_qc_clone_rank_abundance_fig   = vdj_qc_out.qc_figures.flatten().filter { f -> f.name == 'clone_rank_abundance.png'           }.ifEmpty(nofile)
    def vdj_qc_multiple_chains_fig        = vdj_qc_out.qc_figures.flatten().filter { f -> f.name == 'multiple_chains_by_sample.png'      }.ifEmpty(nofile)

    // ── Step 2: pseudobulk → clonotype table (route-specific source) ──────
    def sc_samplesheet = nofile
    if (vdjOnly) {
        VDJ_TO_BULK_SW( vdj_qc_out.contigs_after_qc )
        pseudobulk_map = VDJ_TO_BULK_SW.out.sample_map
        sc_samplesheet = VDJ_TO_BULK_SW.out.samplesheet_utf8
        paired_clones  = VDJ_TO_BULK_SW.out.paired_clones
    } else {
        tcell_out = TCELL_INTEGRATION_SW(
            vdj_qc_out.contigs_after_qc, ch_annotated_obj, ch_project_name
        )
        SC_TO_CDR3_SW( tcell_out.export_cells )
        pseudobulk_map = SC_TO_CDR3_SW.out.sample_map
        paired_clones  = SC_TO_CDR3_SW.out.paired_clones
    }

    // ── Step 2b: merged per-cell TCR object, pre- and post-QC ─────────────
    // Reads the contig tables directly, so it conserves EVERY receptor - including those
    // that never matched a GEX barcode and are therefore absent from the full-SC analysis
    // (which pseudobulks from the post-merge export). Runs after TCELL_INTEGRATION so it
    // can flag which receptors did match.
    mergedvdj_tables = channel.empty()
    if (enabled(params.run_merge_vdj_object)) {
        MERGE_VDJ_OBJECT(
            pickFile(vdj_qc_out.qc_tables, 'contigs_before_qc.tsv', nofile),
            pickFile(vdj_qc_out.qc_tables, 'contigs_after_qc.tsv',  nofile),
            channel.fromPath("${projectDir}/bin/merge_vdj_object.R", checkIfExists: true),
            vdjOnly ? channel.fromPath("${projectDir}/assets/NO_FILE") : tcell_out.export_cells,
            ch_project_name
        )
        mergedvdj_tables = MERGE_VDJ_OBJECT.out.pre_qc.mix(MERGE_VDJ_OBJECT.out.post_qc)
    }

    // ── Step 3: tcrtoolkit pseudobulk QC gate (both routes) ───────────────
    PSEUDOBULK_QC_SW( pseudobulk_map )

    // ── Step 4+5: shared bulk-TCR analysis engine (Option A — no truncation) ──
    // Synthesize the pre-filter-stats sidecar so SC data can use the shared engine's
    // 4-arg SAMPLE step.
    SC_SAMPLE_STATS( PSEUDOBULK_QC_SW.out.sample_map )

    // The report notebooks (QC / discovery brief / details) read a bulk-style samplesheet
    // (sample, alias, patient, timepoint, origin, file). Build one per QC-passing unit from
    // the pseudobulk meta + the user's SC sample sheet. For {sample}__{phenotype} units the
    // base sample is the unit name minus the phenotype suffix.
    def ch_units = PSEUDOBULK_QC_SW.out.sample_map
        .map { meta, f ->
            def unit  = meta.sample.toString()
            def phen  = (meta.phenotype ?: '').toString()
            def base  = (phen && unit.endsWith("__${phen}")) ? unit[0..<(unit.length() - phen.length() - 2)] : unit
            [unit, base, (meta.patient_id ?: '').toString(), phen, f.name].join('\t')
        }
        .collectFile(name: 'report_units.tsv', newLine: true, sort: true,
                     seed: ['unit', 'sample', 'patient_id', 'phenotype', 'file'].join('\t'))
    SC_REPORT_SAMPLESHEET(
        ch_units,
        params.sample_sheet ? file(params.sample_sheet) : nofile,
        vdjOnly ? channel.value(nofile) : tcell_out.export_cells
    )
    def report_samplesheet = SC_REPORT_SAMPLESHEET.out.samplesheet.first()

    // The discovery brief's VDJdb section reads each unit's clonotype table, staged as
    // <sample>_airr.tsv (the same [meta, file] pairs bulk passes from CONVERT).
    def report_clonotype_files = PSEUDOBULK_QC_SW.out.sample_map.collect(flat: false).ifEmpty([])

    BULKTCR_ANALYSIS(
        PSEUDOBULK_QC_SW.out.sample_map,
        SC_SAMPLE_STATS.out.pre_filter_stats,
        PSEUDOBULK_QC_SW.out.concat_cdr3,
        ['sample', 'patient', 'compare'],
        enabled(params.run_reports),
        report_clonotype_files,
        report_samplesheet,
        // Single-cell repertoire metrics count paired clones in cells, not beta CDR3s.
        paired_clones.first(),
        // TCRpheno predicts phenotype from the TCR; with a GEX object it is measured instead.
        vdjOnly
    )

    def concat_cdr3_sorted = BULKTCR_ANALYSIS.out.concat_cdr3_sorted

    // ── Step 6: clonotype export → cluster mapping → consensus → rollups ─────
    // Shared with the bulk modality (CLONOTYPE_CLUSTERING). With a GEX object it maps
    // clusters onto the Seurat; without one it works on a per-cell export synthesised
    // from the pseudobulk.
    conga_report     = channel.empty()
    // Collected tables/ per module for the Master Summary; stay empty where skipped.
    tcell_tables     = channel.empty()
    conga_tables     = channel.empty()
    tcri_tables      = channel.empty()
    tcri_report      = channel.empty()

    if (vdjOnly) {
        BULK_TO_EXPORT( concat_cdr3_sorted, sc_samplesheet, 'cell' )
    } else {
        tcell_tables = tcell_out.tables
    }
    def cluster_export = vdjOnly ? BULK_TO_EXPORT.out.export_cells : tcell_out.export_cells

    // CoNGA needs only the T-cell Seurat (not the sequence clusters), so it runs first: its
    // hit clusters are the fourth method of the consensus on this route.
    def conga_export = channel.value(nofile)
    if (!vdjOnly && enabled(params.run_conga)) {
        conga_out    = CONGA_SW( tcell_out.seurat_tcells_with_tcr, tcell_out.export_cells, ch_project_name )
        conga_report = conga_out.report_html
        conga_tables = conga_out.tables
        conga_export = conga_out.export_cells
    }

    CLONOTYPE_CLUSTERING(
        vdjOnly ? channel.value(nofile) : tcell_out.seurat_tcells_with_tcr,
        !vdjOnly,
        cluster_export,
        BULKTCR_ANALYSIS.out.giana_clusters,
        BULKTCR_ANALYSIS.out.gliph2_cluster_details,
        BULKTCR_ANALYSIS.out.tcrdist_clone_df,
        BULKTCR_ANALYSIS.out.tcrdist_output,
        conga_export,
        ch_project_name
    )
    consensus_report = CLONOTYPE_CLUSTERING.out.consensus_report
    consensus_tables = CLONOTYPE_CLUSTERING.out.consensus_tables
    rep_seurat       = CLONOTYPE_CLUSTERING.out.rep_seurat
    rep_export       = CLONOTYPE_CLUSTERING.out.rep_export

    // ── Step 6b: GEX-only cell-level analyses ─────────────────────────────
    if (!vdjOnly) {
        def enriched_seurat = CLONOTYPE_CLUSTERING.out.enriched_seurat

        // TCRi: immunogenicity scoring on the enriched Seurat. GEX-gated - it needs the
        // transcriptome object, so it cannot run on the VDJ-only route.
        if (enabled(params.run_tcri)) {
            tcri_out    = TCRI_SW( enriched_seurat, tcell_out.export_cells, ch_project_name )
            tcri_tables = tcri_out.tables
            tcri_report = tcri_out.report_html
        }
    }

    // ── Step 7: report notebooks + MASTER_SUMMARY — BOTH routes ───────────
    // Repertoire metrics (diversity, clone-size bins, overlap, sharing) come from the
    // shared engine's REPERTOIRE_METRICS, on paired clones.
    def repertoire_tables = BULKTCR_ANALYSIS.out.repertoire_tables

    // ── Step 7a: report notebooks + their single-cell sections ───────────
    // The QC & Run Info / Discovery Brief / Details notebooks come from the shared engine;
    // each gets its run-level sections and its _sc_*.qmd section, reading the module tables
    // bundled here. The Discovery Brief comes back as markdown for the Master Summary.
    def discovery_md = channel.empty()
    if (enabled(params.run_reports)) {
        REPORT_BUNDLE(
            cluster_export.mix(rep_export).unique().collect(),
            // Per-sample summaries only - the contig-level tables are large and unused.
            vdj_qc_out.qc_tables.flatten()
                .filter { f -> !(f.name.startsWith('contigs_') || f.name == 'clonotypes_all.tsv') }
                .collect().ifEmpty([]),
            channel.empty().mix(PSEUDOBULK_QC_SW.out.qc_summary, PSEUDOBULK_QC_SW.out.v_family)
                .collect().ifEmpty([]),
            channel.empty()
                .mix(BULKTCR_ANALYSIS.out.sample_stats, BULKTCR_ANALYSIS.out.v_family, BULKTCR_ANALYSIS.out.j_family)
                .collect().ifEmpty([]),
            BULKTCR_ANALYSIS.out.shared_cdr3.flatten().collect().ifEmpty([]),
            repertoire_tables.flatten().collect().ifEmpty([]),
            tcell_tables.flatten().collect().ifEmpty([]),
            mergedvdj_tables.flatten().collect().ifEmpty([]),
            consensus_tables.flatten().collect().ifEmpty([]),
            CLONOTYPE_CLUSTERING.out.rollup_tables,
            conga_tables.flatten().collect().ifEmpty([]),
            tcri_tables.flatten().collect().ifEmpty([]),
            vdjOnly ? 'vdj_only' : 'full'
        )
        REPORT( BULKTCR_ANALYSIS.out.reports, report_samplesheet, REPORT_BUNDLE.out.bundle,
                vdjOnly ? 'vdj_only' : 'full' )
        discovery_md = REPORT.out.discovery_md
    }

    if (enabled(params.run_master_summary)) {
        master_barrier = channel.empty()
            .mix(conga_report, consensus_report, tcri_report)
            .collect()
            .ifEmpty([nofile])

        MASTER_SUMMARY_SW(
            rep_seurat,
            rep_export,
            CLONOTYPE_CLUSTERING.out.rollup_tables,
            tcell_tables.flatten().collect().ifEmpty([]),
            conga_tables.flatten().collect().ifEmpty([]),
            consensus_tables.flatten().collect().ifEmpty([]),
            tcri_tables.flatten().collect().ifEmpty([]),
            BULKTCR_ANALYSIS.out.clonal_dynamics_files.map { l -> l ?: [] },
            discovery_md.ifEmpty(nofile).first(),
            master_barrier,
            ch_project_name
        )
    }
}
