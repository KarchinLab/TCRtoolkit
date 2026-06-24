/*
 * SINGLECELL workflow (Scenario 2)
 *
 * Two sub-modes, auto-detected from inputs:
 *
 *   Full SC mode  (--input_annotated_object provided):
 *     VDJ_QC -> TCELL_INTEGRATION
 *       -> Bridge 1 (SC_TO_BULK) -> ANNOTATE -> PATIENT (GIANA + GLIPH2)
 *                                           -> TCRDIST3_MATRIX (per sample)
 *       -> Bridge 2+ (CLUSTER_TO_SC) -> enriched Seurat
 *       -> CONGA -> CONSENSUS -> REPERTOIRE -> MASTER_SUMMARY
 *
 *   VDJ-only mode (--input_annotated_object absent):
 *     VDJ_QC -> Bridge 3 (VDJ_TO_BULK) -> TCRtoolkit bulk analysis
 *     (SC clustering modules are skipped - no GEX object available)
 *
 * Mandatory inputs (both modes):
 *   --input_vdj_contigs   Path glob to Cell Ranger VDJ outs/
 *   --sample_sheet        Sample sheet CSV (columns: sample, path)
 *
 * Additional mandatory input (full SC mode only):
 *   --input_annotated_object   Path to annotated Seurat RDS
 */

include { VDJ_QC_SW }            from '../subworkflows/scratch/vdj_qc.nf'
include { TCELL_INTEGRATION_SW } from '../subworkflows/scratch/tcell_integration.nf'
include { CONGA_SW }             from '../subworkflows/scratch/conga.nf'
include { CONSENSUS_SW }         from '../subworkflows/scratch/consensus_clustering.nf'
include { REPERTOIRE_SW }        from '../subworkflows/scratch/repertoire.nf'
include { MASTER_SUMMARY_SW }    from '../subworkflows/scratch/master_summary.nf'

include { SC_TO_CDR3_SW }        from '../subworkflows/bridges/sc_to_cdr3.nf'
include { VDJ_TO_BULK_SW }       from '../subworkflows/bridges/vdj_to_bulk.nf'
include { CLUSTER_TO_SC_SW }     from '../subworkflows/bridges/cluster_to_sc.nf'

include { ANNOTATE; ANNOTATE_FROM_CONCAT } from '../subworkflows/bulk/annotate.nf'
include { PSEUDOBULK_QC_SW }     from '../subworkflows/bulk/pseudobulk_qc.nf'
include { PATIENT }              from '../subworkflows/bulk/patient.nf'
include { SAMPLE }               from '../subworkflows/bulk/sample.nf'
include { COMPARE }              from '../subworkflows/bulk/compare.nf'

include { TCRDIST3_MATRIX }      from '../modules/bulk/sample/tcrdist3'

workflow SINGLECELL_WORKFLOW {

    def enabled = { x -> x == null || x == true }
    def nofile  = file("${projectDir}/assets/NO_FILE")

    // ── Mandatory inputs (both modes) ─────────────────────────────────────
    if (!params.input_vdj_contigs) error "Please provide --input_vdj_contigs"
    if (!params.sample_sheet)      error "Please provide --sample_sheet"

    def vdjOnly = (!params.input_annotated_object ||
                    params.input_annotated_object == '' ||
                    params.input_annotated_object.endsWith('NO_FILE'))

    if (vdjOnly) {
        log.info "==> VDJ-only mode: no GEX annotated object provided. " +
                 "TCELL_INTEGRATION and SC clustering modules will be skipped."
    }

    ch_sample_sheet  = Channel.fromPath(params.sample_sheet, checkIfExists: true)
    ch_project_name  = Channel.value(params.project_name)
    ch_annotated_obj = vdjOnly
        ? Channel.fromPath("${projectDir}/assets/NO_FILE")
        : Channel.fromPath(params.input_annotated_object, checkIfExists: true)

    // ── Step 1: VDJ QC ────────────────────────────────────────────────────
    vdj_qc_out = VDJ_QC_SW(
        ch_sample_sheet,
        ch_project_name,
        ch_annotated_obj
    )

    // ── Collect VDJ QC tables/figures for Master Summary ─────────────────
    def vdj_qc_per_sample_compact         = vdj_qc_out.qc_tables.flatten().filter { it.name == 'vdj_qc_per_sample_compact.tsv'          }.ifEmpty(nofile)
    def vdj_qc_before_after_summary       = vdj_qc_out.qc_tables.flatten().filter { it.name == 'qc_contigs_before_after_summary.tsv'     }.ifEmpty(nofile)
    def vdj_qc_sample_sheet_resolved      = vdj_qc_out.qc_tables.flatten().filter { it.name == 'sample_sheet_resolved.tsv'               }.ifEmpty(nofile)
    def vdj_qc_clone_rank_abundance       = vdj_qc_out.qc_tables.flatten().filter { it.name == 'clone_rank_abundance.tsv'                }.ifEmpty(nofile)
    def vdj_qc_before_after_retention_fig = vdj_qc_out.qc_figures.flatten().filter { it.name == 'qc_before_after_retention.png'          }.ifEmpty(nofile)
    def vdj_qc_pairing_bar_fig            = vdj_qc_out.qc_figures.flatten().filter { it.name == 'pairing_bar_by_sample.png'              }.ifEmpty(nofile)
    def vdj_qc_clone_rank_abundance_fig   = vdj_qc_out.qc_figures.flatten().filter { it.name == 'clone_rank_abundance.png'               }.ifEmpty(nofile)
    def vdj_qc_multiple_chains_fig        = vdj_qc_out.qc_figures.flatten().filter { it.name == 'multiple_chains_by_sample.png'          }.ifEmpty(nofile)

    // ── Branch: VDJ-only vs Full SC ───────────────────────────────────────
    if (vdjOnly) {

        // ── VDJ-only: Bridge 3 -> TCRtoolkit bulk ──────────────────────────
        VDJ_TO_BULK_SW( vdj_qc_out.contigs_after_qc )
        bulk_sample_map = VDJ_TO_BULK_SW.out.sample_map

        ANNOTATE( bulk_sample_map )

        def levels = params.workflow_level.toLowerCase().tokenize(',')

        if (levels.contains('sample')) {
            SAMPLE( bulk_sample_map, ANNOTATE.out.cdr3_pgen, ANNOTATE.out.olga_stats )
        }
        if (levels.contains('patient')) {
            PATIENT( ANNOTATE.out.processed_samples )
        }
        if (levels.contains('compare')) {
            COMPARE( ANNOTATE.out.concat_cdr3_sorted, ANNOTATE.out.cdr3_pgen )
        }

    } else {

        // ── Full SC mode ───────────────────────────────────────────────────

        // Step 2: T-cell integration
        tcell_out = TCELL_INTEGRATION_SW(
            vdj_qc_out.contigs_after_qc,
            ch_annotated_obj,
            ch_project_name
        )

        // ── Step 3: Bridge 1b - SC export_cells -> CDR3b/TRBV/TRBJ/counts ──
        // Aggregates GEX-annotated per-cell TCR data directly into clonotype
        // counts per sample, skipping the ANNOTATE_PROCESS column-rename step.
        SC_TO_CDR3_SW( tcell_out.export_cells )

        // ── Step 3b: TCRtoolkit QC gate on pseudobulk-derived bulk TCR ────
        PSEUDOBULK_QC_SW( SC_TO_CDR3_SW.out.sample_map )

        // ── Step 4: Sort/dedup/OLGA on QC-passed SC-derived CDR3 data ─────
        ANNOTATE_FROM_CONCAT( PSEUDOBULK_QC_SW.out.sample_map, PSEUDOBULK_QC_SW.out.concat_cdr3 )

        PATIENT( ANNOTATE_FROM_CONCAT.out.processed_samples )

        // TCRdist3 runs per sample directly (no full SAMPLE subworkflow needed)
        TCRDIST3_MATRIX(
            ANNOTATE_FROM_CONCAT.out.processed_samples,
            params.matrix_sparsity,
            params.distance_metric,
            file(params.db_path)
        )

        // ── Step 5: Bridge 2+ - map cluster results -> enriched Seurat ────
        CLUSTER_TO_SC_SW(
            tcell_out.seurat_tcells_with_tcr,
            tcell_out.export_cells,
            PATIENT.out.giana_clusters,
            PATIENT.out.gliph2_cluster_details,
            TCRDIST3_MATRIX.out.clone_df,
            TCRDIST3_MATRIX.out.tcrdist_output.map { _meta, f -> f }
        )

        enriched_seurat = CLUSTER_TO_SC_SW.out.enriched_seurat

        // ── Step 6: CoNGA on enriched Seurat ─────────────────────────────
        if (enabled(params.run_conga)) {
            conga_out = CONGA_SW(
                enriched_seurat,
                tcell_out.export_cells,
                ch_project_name
            )
        }

        // ── Step 7: Consensus clustering ──────────────────────────────────
        if (enabled(params.run_consensus)) {
            def gliph2_export  = CLUSTER_TO_SC_SW.out.gliph2_export.ifEmpty(nofile)
            def tcrdist_export = CLUSTER_TO_SC_SW.out.tcrdist_export.ifEmpty(nofile)
            def giana_export   = CLUSTER_TO_SC_SW.out.giana_export.ifEmpty(nofile)

            CONSENSUS_SW(
                enriched_seurat,
                tcell_out.export_cells,
                gliph2_export,
                tcrdist_export,
                giana_export,
                ch_project_name
            )
        }

        // Determine the best Seurat to pass downstream
        downstream_seurat = enabled(params.run_consensus)
            ? CONSENSUS_SW.out.seurat_with_consensus
            : enriched_seurat

        downstream_export = enabled(params.run_consensus)
            ? CONSENSUS_SW.out.export_cells
            : tcell_out.export_cells

        // ── Step 8: Repertoire ────────────────────────────────────────────
        if (enabled(params.run_repertoire)) {
            REPERTOIRE_SW( downstream_seurat, downstream_export, ch_project_name )
        }

        // ── Step 9: Master summary ────────────────────────────────────────
        if (enabled(params.run_master_summary)) {
            // Barrier: wait for all downstream reports before summary
            done_signals = Channel.empty()
            if (enabled(params.run_conga))      done_signals = done_signals.mix(conga_out.report_html)
            if (enabled(params.run_consensus))  done_signals = done_signals.mix(CONSENSUS_SW.out.report_html)
            if (enabled(params.run_repertoire)) done_signals = done_signals.mix(REPERTOIRE_SW.out.report_html)
            master_barrier = done_signals.collect()

            MASTER_SUMMARY_SW(
                downstream_seurat,
                downstream_export,
                vdj_qc_per_sample_compact,
                vdj_qc_before_after_summary,
                vdj_qc_sample_sheet_resolved,
                vdj_qc_clone_rank_abundance,
                vdj_qc_before_after_retention_fig,
                vdj_qc_pairing_bar_fig,
                vdj_qc_clone_rank_abundance_fig,
                vdj_qc_multiple_chains_fig,
                master_barrier,
                ch_project_name
            )
        }
    }
}
