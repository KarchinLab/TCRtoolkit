/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CLONOTYPE_CLUSTERING

    Shared by both modalities, after the bulk engine's SAMPLE / PATIENT / COMPARE stages:

        CLUSTER_TO_SC  -> map GIANA / GLIPH2 / TCRdist3 clusters onto the export rows
        CONSENSUS      -> reconcile the three methods into one consensus cluster per row
        CLUSTER_ROLLUP -> per-method summary tables for the Master Summary

    `seurat_rds` is the GEX Seurat object on the full single-cell route, and NO_FILE
    otherwise (single-cell VDJ-only route, bulk data). Without it every step runs on the
    export table alone: one row per cell (SC) or per sample x clonotype (bulk).
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CLUSTER_TO_SC_SW } from '../bridges/cluster_to_sc.nf'
include { CONSENSUS_SW }     from '../scratch/consensus_clustering.nf'
include { CLUSTER_ROLLUP }   from '../../modules/bridges/cluster_rollup.nf'

def enabled(x) { x == null || x == true }

workflow CLONOTYPE_CLUSTERING {
    take:
    seurat_rds          // channel: Seurat RDS, or NO_FILE
    has_seurat          // bool
    export_cells        // channel: per-cell / per-clonotype export table
    giana_clusters      // channel: *_giana.txt (per patient)
    gliph2_details      // channel: cluster_member_details (per patient)
    tcrdist_clone_df    // channel: *_clone_df.csv (per sample)
    tcrdist_output      // channel: [meta, distance matrix] (per sample)
    project_name        // value channel

    main:
    def nofile = file("${projectDir}/assets/NO_FILE")

    CLUSTER_TO_SC_SW(
        seurat_rds,
        export_cells,
        giana_clusters,
        gliph2_details,
        tcrdist_clone_df,
        tcrdist_output.map { _meta, f -> f }
    )
    // Without a GEX object CLUSTER_TO_SC writes no Seurat; downstream gets NO_FILE.
    enriched_seurat = has_seurat ? CLUSTER_TO_SC_SW.out.enriched_seurat : channel.value(nofile)

    consensus_report = channel.empty()
    consensus_tables = channel.empty()
    rep_seurat = enriched_seurat
    rep_export = export_cells

    if (enabled(params.run_consensus)) {
        CONSENSUS_SW(
            enriched_seurat,
            export_cells,
            CLUSTER_TO_SC_SW.out.gliph2_export.ifEmpty(nofile),
            CLUSTER_TO_SC_SW.out.tcrdist_export.ifEmpty(nofile),
            CLUSTER_TO_SC_SW.out.giana_export.ifEmpty(nofile),
            project_name
        )
        consensus_report = CONSENSUS_SW.out.report_html
        consensus_tables = CONSENSUS_SW.out.tables
        rep_export       = CONSENSUS_SW.out.export_cells
        if (has_seurat) rep_seurat = CONSENSUS_SW.out.seurat_with_consensus
    }

    // GIANA / GLIPH2 / TCRdist3 write raw per-patient / per-sample output but no rollup
    // table; this summarises them into the files the Master Summary reads.
    CLUSTER_ROLLUP(
        giana_clusters.collect().ifEmpty([]),
        gliph2_details.collect().ifEmpty([]),
        tcrdist_clone_df.collect().ifEmpty([]),
        tcrdist_output.map { _meta, f -> f }.collect().ifEmpty([]),
        rep_export,
        channel.fromPath("${projectDir}/bin/cluster_rollup.py", checkIfExists: true),
        project_name
    )

    rollup_tables = channel.empty()
        .mix(CLUSTER_ROLLUP.out.giana_summary,
             CLUSTER_ROLLUP.out.gliph2_summary,
             CLUSTER_ROLLUP.out.tcrdist3_summary,
             CLUSTER_ROLLUP.out.method_presence,
             CLUSTER_ROLLUP.out.method_cluster_counts,
             CLUSTER_ROLLUP.out.annotation_giana,
             CLUSTER_ROLLUP.out.annotation_gliph2,
             CLUSTER_ROLLUP.out.annotation_tcrdist3,
             CLUSTER_ROLLUP.out.giana_detail,
             CLUSTER_ROLLUP.out.gliph2_detail,
             CLUSTER_ROLLUP.out.tcrdist_detail,
             CLUSTER_ROLLUP.out.giana_vgene,
             CLUSTER_ROLLUP.out.gliph2_vgene)
        .collect().ifEmpty([])

    emit:
    enriched_seurat  = enriched_seurat
    rep_seurat       = rep_seurat
    rep_export       = rep_export
    consensus_report = consensus_report
    consensus_tables = consensus_tables
    rollup_tables    = rollup_tables
}
