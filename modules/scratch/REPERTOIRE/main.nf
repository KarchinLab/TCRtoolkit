process REPERTOIRE {
    tag "${project_name}"
    label 'process_medium'
    container "${params.container}"


    publishDir "${params.outdir}/Repertoire", mode: 'copy', overwrite: true

    input:
    path seurat_rds
    path export_cells
    path qmd
    val  project_name

    output:
    path "Repertoire_Report.html", emit: report_html
    path "Repertoire_Report/tables/*", emit: tables
    path "Repertoire_Report/figures/*", emit: figures

    script:
    """
    mkdir -p Repertoire_Report

    quarto render ${qmd} \
      -P seurat_rds=${seurat_rds} \
      -P export_cells_file=${export_cells} \
      -P outdir=Repertoire_Report \
      -P label_col="${params.label_col}" \
      -P sample_col="${params.sample_col}" \
      -P patient_col="${params.patient_col}" \
      -P condition_col="${params.condition_col}" \
      -P timepoint_col="${params.timepoint_col}" \
      -P batch_col="${params.batch_col}" \
      -P min_cells_per_clone_plot=${params.repertoire_min_cells_per_clone_plot} \
      -P top_n_shared_clones=${params.repertoire_top_n_shared_clones} \
      -P top_n_flux_clones=${params.repertoire_top_n_flux_clones} \
      -P shareability_min_n_groups=${params.repertoire_shareability_min_n_groups} \
      -P min_cells_per_group=${params.repertoire_min_cells_per_group} \
      -P overlap_metric="${params.repertoire_overlap_metric}" \
      -P use_relative_clone_frequencies=${params.repertoire_use_relative_clone_freqs} \
      -P report_label="${project_name} Repertoire"
    """
}