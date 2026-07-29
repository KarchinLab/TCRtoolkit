process CONGA {
    tag "${project_name}"
    label 'process_medium'
    // container "/rsrch8/home/genomic_med/sazaidi/Softwares/SCRATCH_TCR_2025/scratch-tcranalysis2.sif"
    container "${params.container}"

    publishDir "${params.outdir}/CoNGA", mode: 'copy', overwrite: true
//     publishDir "${params.outdir}/CoNGA", mode: params.publish_dir_mode ?: 'copy'    

    input:
      path seurat_rds
      path export_cells
      path qmd
      val  project_name

    output:
      path "CoNGA_Report.html", emit: report_html
      path "CoNGA_Report/data/seurat_with_CoNGA.rds", emit: seurat_with_conga
      path "CoNGA_Report/data/*", emit: data
      path "CoNGA_Report/tables/conga_export_cells.tsv", emit: export_cells
      path "CoNGA_Report/tables/*", emit: tables
      path "CoNGA_Report/figures/*", emit: figures

    script:
      def repo_dir      = params.conga_repo_dir   ?: "/opt/tools/conga"
      def python_bin    = params.conga_python_bin ?: "/opt/conda/envs/tcrenv/bin/python"
      def clone_col     = params.clone_id_col     ?: "clone_id"

    """
    mkdir -p CoNGA_Report

    /opt/quarto/bin/quarto render ${qmd} \
      -P seurat_rds="${seurat_rds}" \
      -P tcr_export_cells_file="${export_cells}" \
      -P filtered_contig_annotations_csvfile="${export_cells}" \
      -P conga_repo_dir="${repo_dir}" \
      -P conga_python_bin="${python_bin}" \
      -P clone_id_col="${clone_col}" \
      -P outdir="CoNGA_Report" \
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
      -P conga_high_cutoff=${params.conga_high_cutoff ?: 0.8} \
      -P conga_mid_cutoff=${params.conga_mid_cutoff ?: 0.5} \
      -P use_quantile_cutoffs_if_score_not_bounded=${params.conga_use_quantile_cutoffs ?: true} \
      -P high_quantile=${params.conga_high_quantile ?: 0.9} \
      -P mid_quantile=${params.conga_mid_quantile ?: 0.5} \
      -P min_cells_per_group=${params.conga_min_cells_per_group ?: 10} \
      -P min_cluster_size_plot=${params.conga_min_cluster_size_plot ?: 5} \
      -P top_n_clusters=${params.conga_top_n_clusters ?: 20} \
      -P max_edges_to_plot=${params.conga_max_edges_to_plot ?: 5000} \
      -P report_label="${project_name} CoNGA"
    """
}

