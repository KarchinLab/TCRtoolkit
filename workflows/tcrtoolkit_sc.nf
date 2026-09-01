/*
 * TCRTOOLKIT_SC — single-cell TCR modality (integrated).
 *
 * Unified spine (see IMPLEMENTATION_SPEC.md §1.3). Both routes converge after their
 * pseudobulk step into: PSEUDOBULK_QC → ANNOTATE → full shared bulk engine
 * (SAMPLE + PATIENT + COMPARE) → repertoire/summary.
 *
 *   Full SC   (--input_annotated_object provided):
 *     VDJ_QC → TCELL_INTEGRATION → SC_TO_CDR3 (phenotype pseudobulk)
 *            → PSEUDOBULK_QC → ANNOTATE
 *            → SAMPLE + PATIENT + COMPARE            (shared engine, full route)
 *            → CLUSTER_TO_SC → CONGA → CONSENSUS     (cell-level, GEX-only)
 *            → REPERTOIRE (cell-level) → MASTER_SUMMARY (full)
 *
 *   VDJ-only  (--input_annotated_object absent):
 *     VDJ_QC → VDJ_TO_BULK (pseudobulk from contigs)
 *            → PSEUDOBULK_QC → ANNOTATE
 *            → SAMPLE + PATIENT + COMPARE            (shared engine, full route)
 *            → BULK_TO_EXPORT (clonotype-level per-cell export)
 *            → REPERTOIRE (clonotype-level) → MASTER_SUMMARY (CoNGA-excluded)
 *            (only CLUSTER_TO_SC / CONGA / CONSENSUS are skipped — they need the GEX substrate)
 */

// ── Single-cell subworkflows (SCRATCH) ────────────────────────────────────
include { VDJ_QC_SW }            from '../subworkflows/scratch/vdj_qc.nf'
include { TCELL_INTEGRATION_SW } from '../subworkflows/scratch/tcell_integration.nf'
include { CONGA_SW }             from '../subworkflows/scratch/conga.nf'
include { CONSENSUS_SW }         from '../subworkflows/scratch/consensus_clustering.nf'
include { REPERTOIRE_SW }        from '../subworkflows/scratch/repertoire.nf'
include { TCRI_SW }             from '../subworkflows/scratch/tcri.nf'
include { MASTER_SUMMARY_SW }    from '../subworkflows/scratch/master_summary.nf'

// ── Bridges (SC ↔ bulk-engine schema conversion) ──────────────────────────
include { SC_TO_CDR3_SW }        from '../subworkflows/bridges/sc_to_cdr3.nf'
include { VDJ_TO_BULK_SW }       from '../subworkflows/bridges/vdj_to_bulk.nf'
include { CLUSTER_TO_SC_SW }     from '../subworkflows/bridges/cluster_to_sc.nf'
include { SC_SAMPLE_STATS }      from '../modules/bridges/sc_sample_stats.nf'
include { BULK_TO_EXPORT }       from '../modules/bridges/bulk_to_export.nf'
include { CLUSTER_ROLLUP }       from '../modules/bridges/cluster_rollup.nf'
include { MERGE_VDJ_OBJECT }     from '../modules/bridges/merge_vdj_object.nf'

// ── Shared bulk-TCR analysis engine (main/local — unchanged behavior) ─────
include { PSEUDOBULK_QC_SW }  from '../subworkflows/local/pseudobulk_qc.nf'
include { BULKTCR_ANALYSIS }  from '../subworkflows/local/bulktcr_analysis.nf'

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
                 "TCELL_INTEGRATION / CLUSTER_TO_SC / CONGA will be skipped."
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

    // ── Step 1b: merged per-cell TCR object, pre- and post-QC ─────────────
    // Real barcodes, all samples pooled. On the VDJ-only route this is the only per-cell
    // object produced - BULK_TO_EXPORT's export is clonotype-level with synthesized cells.
    if (enabled(params.run_merge_vdj_object)) {
        MERGE_VDJ_OBJECT(
            pickFile(vdj_qc_out.qc_tables, 'contigs_before_qc.tsv', nofile),
            pickFile(vdj_qc_out.qc_tables, 'contigs_after_qc.tsv',  nofile),
            channel.fromPath("${projectDir}/bin/merge_vdj_object.R", checkIfExists: true),
            ch_project_name
        )
    }

    // ── Step 2: pseudobulk → clonotype table (route-specific source) ──────
    def sc_samplesheet = nofile
    if (vdjOnly) {
        VDJ_TO_BULK_SW( vdj_qc_out.contigs_after_qc )
        pseudobulk_map = VDJ_TO_BULK_SW.out.sample_map
        sc_samplesheet = VDJ_TO_BULK_SW.out.samplesheet_utf8
    } else {
        tcell_out = TCELL_INTEGRATION_SW(
            vdj_qc_out.contigs_after_qc, ch_annotated_obj, ch_project_name
        )
        SC_TO_CDR3_SW( tcell_out.export_cells )
        pseudobulk_map = SC_TO_CDR3_SW.out.sample_map
    }

    // ── Step 3: tcrtoolkit pseudobulk QC gate (both routes) ───────────────
    PSEUDOBULK_QC_SW( pseudobulk_map )

    // ── Step 4+5: shared bulk-TCR analysis engine (Option A — no truncation) ──
    // Synthesize the pre-filter-stats sidecar so SC data can use the shared engine's
    // 4-arg SAMPLE step. Reports stay off (run_reports: false) — the shared bulk report
    // templates assume bulk-samplesheet metadata (origin/timepoint) that single-cell-derived
    // samplesheets don't reliably carry. Single-cell reporting is handled by REPERTOIRE +
    // MASTER_SUMMARY (below).
    SC_SAMPLE_STATS( PSEUDOBULK_QC_SW.out.sample_map )

    BULKTCR_ANALYSIS(
        PSEUDOBULK_QC_SW.out.sample_map,
        SC_SAMPLE_STATS.out.pre_filter_stats,
        PSEUDOBULK_QC_SW.out.concat_cdr3,
        ['sample', 'patient', 'compare'],
        false,
        channel.value([])
    )

    def concat_cdr3_sorted = BULKTCR_ANALYSIS.out.concat_cdr3_sorted

    // ── Step 6: cell-level clustering — FULL-SC ONLY (needs the GEX/Seurat substrate) ──
    // Produces the enriched/consensus Seurat + export that REPERTOIRE / MASTER_SUMMARY use
    // when a GEX object is present. In VDJ-only mode this whole block is skipped and a
    // clonotype-level export is synthesized instead (below).
    conga_report     = channel.empty()
    consensus_report = channel.empty()
    // Collected tables/ per module for the Master Summary; stay empty where skipped.
    tcell_tables     = channel.empty()
    conga_tables     = channel.empty()
    consensus_tables = channel.empty()
    tcri_tables      = channel.empty()
    tcri_report      = channel.empty()

    if (!vdjOnly) {
        // Reuse tcrdist3 outputs computed inside SAMPLE (no second run).
        CLUSTER_TO_SC_SW(
            tcell_out.seurat_tcells_with_tcr,
            tcell_out.export_cells,
            BULKTCR_ANALYSIS.out.giana_clusters,
            BULKTCR_ANALYSIS.out.gliph2_cluster_details,
            BULKTCR_ANALYSIS.out.tcrdist_clone_df,
            BULKTCR_ANALYSIS.out.tcrdist_output.map { _meta, f -> f }
        )
        enriched_seurat = CLUSTER_TO_SC_SW.out.enriched_seurat
        tcell_tables    = tcell_out.tables

        if (enabled(params.run_conga)) {
            conga_out    = CONGA_SW( enriched_seurat, tcell_out.export_cells, ch_project_name )
            conga_report = conga_out.report_html
            conga_tables = conga_out.tables
        }

        // TCRi: immunogenicity scoring on the enriched Seurat. GEX-gated - it needs the
        // transcriptome object, so it cannot run on the VDJ-only route.
        if (enabled(params.run_tcri)) {
            tcri_out    = TCRI_SW( enriched_seurat, tcell_out.export_cells, ch_project_name )
            tcri_tables = tcri_out.tables
            tcri_report = tcri_out.report_html
        }

        if (enabled(params.run_consensus)) {
            def gliph2_export  = CLUSTER_TO_SC_SW.out.gliph2_export.ifEmpty(nofile)
            def tcrdist_export = CLUSTER_TO_SC_SW.out.tcrdist_export.ifEmpty(nofile)
            def giana_export   = CLUSTER_TO_SC_SW.out.giana_export.ifEmpty(nofile)
            CONSENSUS_SW(
                enriched_seurat, tcell_out.export_cells,
                gliph2_export, tcrdist_export, giana_export, ch_project_name
            )
            consensus_report = CONSENSUS_SW.out.report_html
            consensus_tables = CONSENSUS_SW.out.tables
        }

        rep_seurat = enabled(params.run_consensus) ? CONSENSUS_SW.out.seurat_with_consensus : enriched_seurat
        rep_export = enabled(params.run_consensus) ? CONSENSUS_SW.out.export_cells          : tcell_out.export_cells

    } else {
        // VDJ-only: no Seurat. Build a clonotype-level per-cell export from the pseudobulk
        // so REPERTOIRE / MASTER_SUMMARY run without a GEX object (CoNGA/consensus skipped).
        BULK_TO_EXPORT( concat_cdr3_sorted, sc_samplesheet )
        rep_seurat = channel.fromPath("${projectDir}/assets/NO_FILE")
        rep_export = BULK_TO_EXPORT.out.export_cells
    }

    // ── Step 7: REPERTOIRE + MASTER_SUMMARY — BOTH routes ─────────────────
    // Cell-level + full with a GEX object; clonotype-level repertoire and a CoNGA-excluded
    // summary without one (only CoNGA and the cell-cluster mapping are truly GEX-gated).
    repertoire_report = channel.empty()
    repertoire_tables = channel.empty()
    if (enabled(params.run_repertoire)) {
        REPERTOIRE_SW( rep_seurat, rep_export, ch_project_name )
        repertoire_report = REPERTOIRE_SW.out.report_html
        repertoire_tables = REPERTOIRE_SW.out.tables
    }

    // ── Step 7b: rollups for GIANA / GLIPH2 / TCRdist3 ────────────────────
    // These three write raw per-patient / per-sample output but no rollup table, so the
    // Master Summary had nothing to read for them and they always reported as absent.
    CLUSTER_ROLLUP(
        BULKTCR_ANALYSIS.out.giana_clusters.collect().ifEmpty([]),
        BULKTCR_ANALYSIS.out.gliph2_cluster_details.collect().ifEmpty([]),
        BULKTCR_ANALYSIS.out.tcrdist_clone_df.collect().ifEmpty([]),
        BULKTCR_ANALYSIS.out.tcrdist_output.map { _meta, f -> f }.collect().ifEmpty([]),
        rep_export,
        channel.fromPath("${projectDir}/bin/cluster_rollup.py", checkIfExists: true),
        ch_project_name
    )

    if (enabled(params.run_master_summary)) {
        master_barrier = channel.empty()
            .mix(conga_report, consensus_report, repertoire_report, tcri_report)
            .collect()
            .ifEmpty([nofile])

        // Aggregated per-sample stats, built here rather than emitted from the shared
        // engine so no bulk-side file needs changing.
        def sample_stats_agg = BULKTCR_ANALYSIS.out.sample_csv
            .collectFile(name: "sample_stats.csv", keepHeader: true, skip: 1, sort: true)

        def rollup_tables = channel.empty()
            .mix(CLUSTER_ROLLUP.out.giana_summary,
                 CLUSTER_ROLLUP.out.gliph2_summary,
                 CLUSTER_ROLLUP.out.tcrdist3_summary,
                 CLUSTER_ROLLUP.out.method_presence,
                 CLUSTER_ROLLUP.out.method_cluster_counts,
                 CLUSTER_ROLLUP.out.annotation_giana,
                 CLUSTER_ROLLUP.out.annotation_gliph2,
                 CLUSTER_ROLLUP.out.annotation_tcrdist3)
            .collect().ifEmpty([])

        def pseudobulk_tables = channel.empty()
            .mix(PSEUDOBULK_QC_SW.out.qc_summary, PSEUDOBULK_QC_SW.out.v_family)
            .collect().ifEmpty([])

        def sample_tables = channel.empty()
            .mix(sample_stats_agg, BULKTCR_ANALYSIS.out.v_family, BULKTCR_ANALYSIS.out.j_family)
            .collect().ifEmpty([])

        MASTER_SUMMARY_SW(
            rep_seurat,
            rep_export,
            vdj_qc_out.qc_tables.flatten().collect().ifEmpty([]),
            pseudobulk_tables,
            sample_tables,
            BULKTCR_ANALYSIS.out.shared_cdr3.flatten().collect().ifEmpty([]),
            rollup_tables,
            repertoire_tables.flatten().collect().ifEmpty([]),
            tcell_tables.flatten().collect().ifEmpty([]),
            conga_tables.flatten().collect().ifEmpty([]),
            consensus_tables.flatten().collect().ifEmpty([]),
            tcri_tables.flatten().collect().ifEmpty([]),
            master_barrier,
            ch_project_name
        )
    }
}
