process TCRI {

    tag "Running TCRi - ${project_name}"
    label 'process_medium'
    container "${params.container}"

    publishDir "${params.outdir}/TCRi", mode: 'copy', overwrite: true

    input:
        path qmd
        path seurat_rds
        path export_cells
        val  project_name

    output:
        path "TCRi_Report.html", emit: report_html
        path "TCRi_Report/data/seurat_with_TCRi.rds", emit: seurat_with_tcri
        path "TCRi_Report/data/*", emit: data
        path "TCRi_Report/tables/tcri_export_cells.tsv", emit: export_cells
        path "TCRi_Report/tables/*", emit: tables
        path "TCRi_Report/figures/*", emit: figures

    script:
        def tcri_scores = params.tcri_scores_file ?: ""

    """
    mkdir -p TCRi_Report
    mkdir -p .cache/quarto
    export XDG_CACHE_HOME="\$PWD/.cache"
    export QUARTO_CACHE_DIR="\$PWD/.cache/quarto"
    export XDG_DATA_HOME="\$PWD/.cache"
    export QUARTO_PRINT_STACK=true
    export HOME="\$PWD"

    export LD_LIBRARY_PATH="/opt/conda/envs/tcrenv/lib:\$LD_LIBRARY_PATH"
    export RETICULATE_PYTHON="/opt/conda/envs/tcrenv/bin/python"

    echo "Testing Python environment natively..."
    /opt/conda/envs/tcrenv/bin/python -c "import tcri; print('PYTHON IMPORT SUCCESS!')"
    echo "Starting Quarto Render..."

    quarto render ${qmd} \
      -P python_bin="/opt/conda/envs/tcrenv/bin/python" \
      -P seurat_rds="${seurat_rds}" \
      -P tcr_export_cells_file="${export_cells}" \
      -P tcri_scores_file="${tcri_scores}" \
      -P outdir="TCRi_Report" \
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
      -P tcri_high_cutoff=${params.tcri_high_cutoff} \
      -P tcri_mid_cutoff=${params.tcri_mid_cutoff} \
      -P use_quantile_cutoffs_if_score_not_bounded=${params.tcri_use_quantile_cutoffs} \
      -P high_quantile=${params.tcri_high_quantile} \
      -P mid_quantile=${params.tcri_mid_quantile} \
      -P min_cells_per_group=${params.tcri_min_cells_per_group} \
      -P min_clone_size_for_assoc=${params.tcri_min_clone_size_for_assoc} \
      -P top_n_high_tcri_clones=${params.tcri_top_n_high_clones} \
      -P report_label="${project_name} TCRi"
    """
}

