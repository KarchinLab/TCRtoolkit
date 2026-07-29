/*
 * Bridge 2+: CLUSTER_TO_SC
 *
 * Maps bulk GIANA / GLIPH2 / TCRdist3 cluster assignments back onto single
 * cells via CDR3b as the join key.
 *
 * All cluster files are discovered by pattern in the working directory:
 *   *_giana.txt                    - GIANA per-patient output
 *   cluster_member_details*.txt    - GLIPH2 per-patient output
 *   *_clone_df.csv                 - TCRdist3 per-sample clone table
 *   *_distance_matrix.*            - TCRdist3 per-sample distance matrix
 *
 * Outputs:
 *   enriched_seurat.rds        - Seurat with giana/gliph2/tcrdist cluster columns
 *   giana_export_cells.tsv     - per-cell TSV with giana_cluster column
 *   gliph2_export_cells.tsv    - per-cell TSV with gliph2_cluster column
 *   tcrdist_export_cells.tsv   - per-cell TSV with tcrdist_cluster column
 */

process CLUSTER_TO_SC {
    tag "Bridge 2+ - Cluster to SC annotation"
    label 'process_medium'
    container "${params.container}"

    publishDir "${params.outdir}/bridge/cluster_to_sc", mode: 'copy', overwrite: true

    input:
    path seurat_rds
    path export_cells
    path giana_files        // collected *_giana.txt files (one per patient)
    path gliph2_files       // collected cluster_member_details dirs/files
    path tcrdist_clone_dfs  // collected *_clone_df.csv files (one per sample)
    path tcrdist_matrices   // collected *_distance_matrix.* files (one per sample)
    val  tcrdist_radius

    output:
    path "enriched_seurat.rds",       emit: enriched_seurat
    path "giana_export_cells.tsv",    emit: giana_export,    optional: true
    path "gliph2_export_cells.tsv",   emit: gliph2_export,   optional: true
    path "tcrdist_export_cells.tsv",  emit: tcrdist_export,  optional: true

    script:
    """
    Rscript ${projectDir}/bin/enrich_seurat.R \\
        --seurat_rds     "${seurat_rds}" \\
        --export_cells   "${export_cells}" \\
        --tcrdist_radius ${tcrdist_radius}
    """
}
