#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { CLUSTER_TO_SC } from '../../modules/bridges/cluster_to_sc.nf'

/*
 * CLUSTER_TO_SC_SW  (Bridge 2+)
 *
 * Collects per-patient GIANA / GLIPH2 and per-sample TCRdist3 outputs, then
 * maps all cluster assignments back onto single cells via the Seurat object.
 *
 * Takes:
 *   seurat_rds        - seurat_tcells_with_tcr from TCELL_INTEGRATION_SW
 *   export_cells      - export_cells.tsv from TCELL_INTEGRATION_SW
 *   giana_clusters    - channel of *_giana.txt files from PATIENT.out.giana_clusters
 *   gliph2_details    - channel of cluster_member_details files from
 *                       PATIENT.out.gliph2_cluster_details  (may be empty)
 *   tcrdist_clone_dfs - channel of *_clone_df.csv from TCRDIST3_MATRIX.out.clone_df
 *   tcrdist_matrices  - channel of distance matrix files from
 *                       TCRDIST3_MATRIX.out.tcrdist_output (file only, meta stripped)
 *
 * Emits:
 *   enriched_seurat   - Seurat RDS with giana/gliph2/tcrdist cluster columns added
 *   giana_export      - per-cell TSV with giana_cluster column
 *   gliph2_export     - per-cell TSV with gliph2_cluster column
 *   tcrdist_export    - per-cell TSV with tcrdist_cluster column
 */
workflow CLUSTER_TO_SC_SW {
    take:
    seurat_rds
    export_cells
    giana_clusters
    gliph2_details
    tcrdist_clone_dfs
    tcrdist_matrices

    main:
    def nofile = file("${projectDir}/assets/NO_FILE")

    // Collect all per-patient / per-sample files into single inputs
    giana_collected    = giana_clusters.collect().ifEmpty([nofile])
    gliph2_collected   = gliph2_details.collect().ifEmpty([nofile])
    clone_dfs_collected = tcrdist_clone_dfs.collect().ifEmpty([nofile])
    matrices_collected  = tcrdist_matrices.collect().ifEmpty([nofile])

    CLUSTER_TO_SC(
        seurat_rds,
        export_cells,
        giana_collected,
        gliph2_collected,
        clone_dfs_collected,
        matrices_collected,
        params.tcrdist_radius ?: 24
    )

    emit:
    enriched_seurat = CLUSTER_TO_SC.out.enriched_seurat
    giana_export    = CLUSTER_TO_SC.out.giana_export
    gliph2_export   = CLUSTER_TO_SC.out.gliph2_export
    tcrdist_export  = CLUSTER_TO_SC.out.tcrdist_export
}
