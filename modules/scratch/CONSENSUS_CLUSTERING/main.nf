process CONSENSUS_CLUSTERING {
    tag "${project_name}"
    label 'process_medium'
    container "${params.container}"

    publishDir "${params.outdir}/Consensus_Clustering", mode: 'copy', overwrite: true

    input:
      path seurat_rds
      path export_cells
      path gliph_export_cells
      path tcrdist_export_cells
      path giana_export_cells
      path qmd
      val  project_name

    output:
      path "Clonotype_Clustering_Consensus_Report.html", emit: report_html
      path "Clonotype_Clustering_Consensus_Report/data/seurat_with_consensus_clonotype_clusters.rds", emit: seurat_with_consensus
      path "Clonotype_Clustering_Consensus_Report/tables/consensus_export_cells.tsv", emit: export_cells
      path "Clonotype_Clustering_Consensus_Report/tables/*", emit: tables
      path "Clonotype_Clustering_Consensus_Report/figures/*", emit: figures

    script:
    """
    mkdir -p Clonotype_Clustering_Consensus_Report

    quarto render ${qmd} \\
      -P seurat_rds=${seurat_rds} \\
      -P export_cells_file=${export_cells} \\
      -P gliph_export_cells_file=${gliph_export_cells} \\
      -P tcrdist_export_cells_file=${tcrdist_export_cells} \\
      -P giana_export_cells_file=${giana_export_cells} \\
      -P outdir=Clonotype_Clustering_Consensus_Report \\
      -P label_col="${params.label_col}" \\
      -P sample_col="${params.sample_col}" \\
      -P patient_col="${params.patient_col}" \\
      -P condition_col="${params.condition_col}" \\
      -P timepoint_col="${params.timepoint_col}" \\
      -P batch_col="${params.batch_col}" \\
      -P reduction_use="${params.reduction_use}" \\
      -P make_umap_if_missing=${params.make_umap_if_missing} \\
      -P umap_dims_max=${params.umap_dims_max} \\
      -P umap_nfeatures=${params.umap_nfeatures} \\
      -P raster_large_umap=${params.raster_large_umap} \\
      -P consensus_min_methods=${params.consensus_min_methods} \\
      -P use_majority_vote=${params.consensus_use_majority_vote} \\
      -P assign_singleton_consensus=${params.consensus_assign_singleton} \\
      -P consensus_label_prefix="${params.consensus_label_prefix}" \\
      -P min_cells_per_group=${params.consensus_min_cells_per_group} \\
      -P min_cluster_size_plot=${params.consensus_min_cluster_size_plot} \\
      -P top_n_clusters=${params.consensus_top_n_clusters} \\
      -P report_label="${project_name} Consensus Clustering"
    """
}
