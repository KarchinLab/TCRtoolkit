process TCELL_INTEGRATION {

    tag "${project_name}"
    label 'process_medium'

    container "${params.container}"

    publishDir "${params.outdir}/TCell_Integration",mode: 'copy', overwrite: true
    // publishDir "${params.outdir}/${params.project_name}", 
    // mode: 'copy', overwrite: true


    input:
      path (contigs_after_qc)
      path (annotated_object)
      path (qmd)
      val  (project_name)

    output:
      path "TCell_Integration_Report.html", emit: report_html
      path "TCell_Integration_Report/data/seurat_tcells_with_TCR.rds", emit: seurat_tcells_with_tcr
      path "TCell_Integration_Report/tables/tcr_export_cells_with_embedding.tsv", emit: export_cells
      path("TCell_Integration_Report/tables/*"),                    emit: tables
      path("TCell_Integration_Report/figures/*"),                   emit: figures




    script:
    """
    mkdir -p TCell_Integration_Report

    quarto render ${qmd} \
      -P contigs_file=${contigs_after_qc} \
      -P seurat_rds=${annotated_object} \
      -P outdir=TCell_Integration_Report \
      -P label_col="${params.label_col}" \
      -P sample_col="${params.sample_col}" \
      -P patient_col="${params.patient_col}" \
      -P condition_col="${params.condition_col}" \
      -P timepoint_col="${params.timepoint_col}" \
      -P batch_col="${params.batch_col}" \
      -P cells_mode="${params.cells_mode}" \
      -P filter_to_t_ab=${params.filter_to_t_ab} \
      -P clone_call_preference="${params.clone_call_preference}" \
      -P keep_na_clonotypes=${params.keep_na_clonotypes} \
      -P harmonize_apply_to="${params.harmonize_apply_to}" \
      -P harmonization_overlap_threshold=${params.harmonization_overlap_threshold} \
      -P minimum_final_overlap_fraction=${params.minimum_final_overlap_fraction} \
      -P subset_tcells=${params.subset_tcells} \
      -P tcell_regex="${params.tcell_regex}" \
      -P reduction_use="${params.reduction_use}" \
      -P make_umap_if_missing=${params.make_umap_if_missing} \
      -P umap_dims_max=${params.umap_dims_max} \
      -P umap_nfeatures=${params.umap_nfeatures} \
      -P raster_large_umap=${params.raster_large_umap} \
      -P min_clone_size_plot=${params.min_clone_size_plot} \
      -P top_n_clones_umap=${params.top_n_clones_umap} \
      -P top_n_clone_table=${params.top_n_clone_table} \
      -P report_label="${project_name} TCell Integration"
    """
}