process GIANA {
    tag "${project_name}"
    label 'process_medium'
    container "${params.container}"
    
    publishDir "${params.outdir}/GIANA", mode: 'copy', overwrite: true

    input:
    path seurat_rds
    path export_cells
    path qmd
    val  project_name

    output:
    path "GIANA_Report.html", emit: report_html
    path "GIANA_Report/data/seurat_with_GIANA.rds", emit: seurat_with_giana
    path "GIANA_Report/tables/giana_export_cells.tsv", emit: export_cells
    path "GIANA_Report/tables/*", emit: tables
    path "GIANA_Report/figures/*", emit: figures

    script:
    def giana_clusters  = params.giana_clusters_file ?: ""
    def giana_summary   = params.giana_summary_file ?: ""
    def giana_neighbors = params.giana_neighbors_file ?: ""

    """
    mkdir -p GIANA_Report

    quarto render ${qmd} \
      -P seurat_rds="${seurat_rds}" \
      -P tcr_export_cells_file="${export_cells}" \
      -P giana_clusters_file="${giana_clusters}" \
      -P giana_summary_file="${giana_summary}" \
      -P giana_neighbors_file="${giana_neighbors}" \
      -P outdir="GIANA_Report" \
      -P label_col="${params.label_col}" \
      -P sample_col="${params.sample_col}" \
      -P patient_col="${params.patient_col}" \
      -P condition_col="${params.condition_col}" \
      -P timepoint_col="${params.timepoint_col}" \
      -P batch_col="${params.batch_col}" \
      -P reduction_use="${params.reduction_use}" \
      -P make_umap_if_missing=${params.make_umap_if_missing} \
      -P umap_dims_max=${params.umap_dims_max} \
      -P umap_nfeatures=${params.umap_nfeatures} \
      -P raster_large_umap=${params.raster_large_umap} \
      -P clone_match_mode="${params.giana_clone_match_mode}" \
      -P derive_clone_key_from_CTaa=${params.giana_derive_clone_key_from_CTaa} \
      -P min_cells_per_group=${params.giana_min_cells_per_group} \
      -P min_cluster_size_plot=${params.giana_min_cluster_size_plot} \
      -P top_n_clusters=${params.giana_top_n_clusters} \
      -P report_label="${project_name} GIANA"
    """
}
