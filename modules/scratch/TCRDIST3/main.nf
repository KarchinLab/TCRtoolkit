process TCRDIST3 {
    tag "${project_name}"
    label 'process_medium'
    container "${params.container}"

    

    publishDir "${params.outdir}/TCRdist3", mode: 'copy', overwrite: true

    input:
    path seurat_rds
    path export_cells
    path qmd
    val  project_name

    output:
    path "TCRdist3_Report.html", emit: report_html
    path "TCRdist3_Report/data/seurat_with_TCRdist3.rds", emit: seurat_with_tcrdist3
    path "TCRdist3_Report/tables/tcrdist3_export_cells.tsv", emit: export_cells
    path "TCRdist3_Report/tables/*", emit: tables
    path "TCRdist3_Report/figures/*", emit: figures

    script:
    def tcrdist_clusters  = params.tcrdist_clusters_file ?: ""
    def tcrdist_neighbors = params.tcrdist_neighbors_file ?: ""
    def tcrdist_distance  = params.tcrdist_distance_file ?: ""

    """
    mkdir -p TCRdist3_Report
    mkdir -p .cache/quarto
    export XDG_CACHE_HOME="\$PWD/.cache"
    export QUARTO_CACHE_DIR="\$PWD/.cache/quarto"
    export XDG_DATA_HOME="\$PWD/.cache"
    export QUARTO_PRINT_STACK=true
    export HOME="\$PWD"

    export LD_LIBRARY_PATH="/opt/conda/envs/tcrenv/lib:\$LD_LIBRARY_PATH"
    export RETICULATE_PYTHON="/opt/conda/envs/tcrenv/bin/python"

    echo "Testing Python environment natively..."
    /opt/conda/envs/tcrenv/bin/python -c "import llvmlite; import numba; print('PYTHON IMPORT SUCCESS!')"
    echo "Starting Quarto Render..."


    quarto render ${qmd} \
      -P seurat_rds="${seurat_rds}" \
      -P tcr_export_cells_file="${export_cells}" \
      -P tcrdist_clusters_file="${tcrdist_clusters}" \
      -P tcrdist_neighbors_file="${tcrdist_neighbors}" \
      -P tcrdist_distance_file="${tcrdist_distance}" \
      -P outdir="TCRdist3_Report" \
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
      -P clone_match_mode="${params.tcrdist_clone_match_mode}" \
      -P derive_clone_key_from_CTaa=${params.tcrdist_derive_clone_key_from_CTaa} \
      -P min_cells_per_group=${params.tcrdist_min_cells_per_group} \
      -P min_cluster_size_plot=${params.tcrdist_min_cluster_size_plot} \
      -P top_n_clusters=${params.tcrdist_top_n_clusters} \
      -P top_n_neighbors=${params.tcrdist_top_n_neighbors} \
      -P report_label="${project_name} TCRdist3"
    """
}
