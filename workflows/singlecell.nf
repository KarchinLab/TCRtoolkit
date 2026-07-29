/*
 * SINGLECELL_WORKFLOW — single-cell TCR modality (integrated).
 *
 * Unified spine (see IMPLEMENTATION_SPEC.md §1.3). Both routes converge after their
 * pseudobulk step into: PSEUDOBULK_QC → ANNOTATE_FROM_CONCAT → full shared bulk engine
 * (SAMPLE + PATIENT + COMPARE) → repertoire/summary.
 *
 *   Full SC   (--input_annotated_object provided):
 *     VDJ_QC → TCELL_INTEGRATION → SC_TO_CDR3 (phenotype pseudobulk)
 *            → PSEUDOBULK_QC → ANNOTATE_FROM_CONCAT
 *            → SAMPLE + PATIENT + COMPARE            (shared engine, full route)
 *            → CLUSTER_TO_SC → CONGA → CONSENSUS     (cell-level, GEX-only)
 *            → REPERTOIRE (cell-level) → MASTER_SUMMARY (full)
 *
 *   VDJ-only  (--input_annotated_object absent):
 *     VDJ_QC → VDJ_TO_BULK (pseudobulk from contigs)
 *            → PSEUDOBULK_QC → ANNOTATE_FROM_CONCAT
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
include { MASTER_SUMMARY_SW }    from '../subworkflows/scratch/master_summary.nf'

// ── Bridges (SC ↔ bulk-engine schema conversion) ──────────────────────────
include { SC_TO_CDR3_SW }        from '../subworkflows/bridges/sc_to_cdr3.nf'
include { VDJ_TO_BULK_SW }       from '../subworkflows/bridges/vdj_to_bulk.nf'
include { CLUSTER_TO_SC_SW }     from '../subworkflows/bridges/cluster_to_sc.nf'
include { SC_SAMPLE_STATS }      from '../modules/bridges/sc_sample_stats.nf'
include { BULK_TO_EXPORT }       from '../modules/bridges/bulk_to_export.nf'

// ── Shared bulk engine (main/local — unchanged behavior) ──────────────────
include { ANNOTATE_FROM_CONCAT } from '../subworkflows/local/annotate.nf'
include { PSEUDOBULK_QC_SW }     from '../subworkflows/local/pseudobulk_qc.nf'
include { SAMPLE }               from '../subworkflows/local/sample.nf'
include { PATIENT }              from '../subworkflows/local/patient.nf'
include { COMPARE }              from '../subworkflows/local/compare.nf'


workflow SINGLECELL_WORKFLOW {

    def enabled = { x -> x == null || x == true }
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

    ch_sample_sheet  = Channel.fromPath(params.sample_sheet, checkIfExists: true)
    ch_project_name  = Channel.value(params.project_name)
    ch_annotated_obj = vdjOnly
        ? Channel.fromPath("${projectDir}/assets/NO_FILE")
        : Channel.fromPath(params.input_annotated_object, checkIfExists: true)

    // ── Step 1: VDJ QC (both routes) ──────────────────────────────────────
    vdj_qc_out = VDJ_QC_SW( ch_sample_sheet, ch_project_name, ch_annotated_obj )

    // VDJ QC tables/figures kept for the Master Summary (full-SC route)
    def vdj_qc_per_sample_compact         = vdj_qc_out.qc_tables.flatten().filter { it.name == 'vdj_qc_per_sample_compact.tsv'      }.ifEmpty(nofile)
    def vdj_qc_before_after_summary       = vdj_qc_out.qc_tables.flatten().filter { it.name == 'qc_contigs_before_after_summary.tsv' }.ifEmpty(nofile)
    def vdj_qc_sample_sheet_resolved      = vdj_qc_out.qc_tables.flatten().filter { it.name == 'sample_sheet_resolved.tsv'           }.ifEmpty(nofile)
    def vdj_qc_clone_rank_abundance       = vdj_qc_out.qc_tables.flatten().filter { it.name == 'clone_rank_abundance.tsv'            }.ifEmpty(nofile)
    def vdj_qc_before_after_retention_fig = vdj_qc_out.qc_figures.flatten().filter { it.name == 'qc_before_after_retention.png'      }.ifEmpty(nofile)
    def vdj_qc_pairing_bar_fig            = vdj_qc_out.qc_figures.flatten().filter { it.name == 'pairing_bar_by_sample.png'          }.ifEmpty(nofile)
    def vdj_qc_clone_rank_abundance_fig   = vdj_qc_out.qc_figures.flatten().filter { it.name == 'clone_rank_abundance.png'           }.ifEmpty(nofile)
    def vdj_qc_multiple_chains_fig        = vdj_qc_out.qc_figures.flatten().filter { it.name == 'multiple_chains_by_sample.png'      }.ifEmpty(nofile)

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

    // ── Step 4: OLGA annotation on QC-passed pseudobulk (shared engine) ───
    ANNOTATE_FROM_CONCAT( PSEUDOBULK_QC_SW.out.sample_map, PSEUDOBULK_QC_SW.out.concat_cdr3 )

    def processed_samples  = ANNOTATE_FROM_CONCAT.out.processed_samples
    def concat_cdr3_sorted = ANNOTATE_FROM_CONCAT.out.concat_cdr3_sorted
    def cdr3_pgen          = ANNOTATE_FROM_CONCAT.out.cdr3_pgen
    def olga_stats         = ANNOTATE_FROM_CONCAT.out.olga_stats

    // ── Step 5: FULL shared bulk route (Option A — no truncation) ─────────
    // Synthesize the pre-filter-stats sidecar so SC data can use main's 4-arg SAMPLE.
    SC_SAMPLE_STATS( processed_samples )

    SAMPLE( processed_samples, SC_SAMPLE_STATS.out.pre_filter_stats, cdr3_pgen, olga_stats )
    PATIENT( processed_samples )
    COMPARE( concat_cdr3_sorted, cdr3_pgen )

    // Note: the shared bulk template_qc REPORT is intentionally NOT run here — it renders a
    // bulk-samplesheet-metadata QC notebook that doesn't fit single-cell-derived data.
    // Single-cell reporting is handled by REPERTOIRE + MASTER_SUMMARY (below).

    // ── Step 6: cell-level clustering — FULL-SC ONLY (needs the GEX/Seurat substrate) ──
    // Produces the enriched/consensus Seurat + export that REPERTOIRE / MASTER_SUMMARY use
    // when a GEX object is present. In VDJ-only mode this whole block is skipped and a
    // clonotype-level export is synthesized instead (below).
    conga_report     = Channel.empty()
    consensus_report = Channel.empty()

    if (!vdjOnly) {
        // Reuse tcrdist3 outputs computed inside SAMPLE (no second run).
        CLUSTER_TO_SC_SW(
            tcell_out.seurat_tcells_with_tcr,
            tcell_out.export_cells,
            PATIENT.out.giana_clusters,
            PATIENT.out.gliph2_cluster_details,
            SAMPLE.out.tcrdist_clone_df,
            SAMPLE.out.tcrdist_output.map { _meta, f -> f }
        )
        enriched_seurat = CLUSTER_TO_SC_SW.out.enriched_seurat

        if (enabled(params.run_conga)) {
            conga_out    = CONGA_SW( enriched_seurat, tcell_out.export_cells, ch_project_name )
            conga_report = conga_out.report_html
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
        }

        rep_seurat = enabled(params.run_consensus) ? CONSENSUS_SW.out.seurat_with_consensus : enriched_seurat
        rep_export = enabled(params.run_consensus) ? CONSENSUS_SW.out.export_cells          : tcell_out.export_cells

    } else {
        // VDJ-only: no Seurat. Build a clonotype-level per-cell export from the pseudobulk
        // so REPERTOIRE / MASTER_SUMMARY run without a GEX object (CoNGA/consensus skipped).
        BULK_TO_EXPORT( concat_cdr3_sorted, sc_samplesheet )
        rep_seurat = Channel.fromPath("${projectDir}/assets/NO_FILE")
        rep_export = BULK_TO_EXPORT.out.export_cells
    }

    // ── Step 7: REPERTOIRE + MASTER_SUMMARY — BOTH routes ─────────────────
    // Cell-level + full with a GEX object; clonotype-level repertoire and a CoNGA-excluded
    // summary without one (only CoNGA and the cell-cluster mapping are truly GEX-gated).
    repertoire_report = Channel.empty()
    if (enabled(params.run_repertoire)) {
        REPERTOIRE_SW( rep_seurat, rep_export, ch_project_name )
        repertoire_report = REPERTOIRE_SW.out.report_html
    }

    if (enabled(params.run_master_summary)) {
        master_barrier = Channel.empty()
            .mix(conga_report, consensus_report, repertoire_report)
            .collect()
            .ifEmpty([nofile])

        MASTER_SUMMARY_SW(
            rep_seurat, rep_export,
            vdj_qc_per_sample_compact, vdj_qc_before_after_summary,
            vdj_qc_sample_sheet_resolved, vdj_qc_clone_rank_abundance,
            vdj_qc_before_after_retention_fig, vdj_qc_pairing_bar_fig,
            vdj_qc_clone_rank_abundance_fig, vdj_qc_multiple_chains_fig,
            master_barrier, ch_project_name
        )
    }
}
