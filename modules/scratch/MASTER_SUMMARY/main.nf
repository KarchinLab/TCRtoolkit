process MASTER_SUMMARY {
    tag "${project_name}"
    label 'process_medium'
    container "${params.container}"

    publishDir "${params.outdir}/Master_Summary", mode: 'copy', overwrite: true

    input:
      path seurat_rds
      path export_cells

      path vdj_qc_per_sample_compact
      path vdj_qc_before_after_summary
      path vdj_qc_sample_sheet_resolved
      path vdj_qc_clone_rank_abundance

      path vdj_qc_before_after_retention_fig
      path vdj_qc_pairing_bar_fig
      path vdj_qc_clone_rank_abundance_fig
      path vdj_qc_multiple_chains_fig

      path qmd
      val  barrier_done
      val  project_name

    output:
      path "Master_Summary_Report.html", emit: report_html
      path "Master_Summary_Report/tables/*", emit: tables, optional: true
      path "Master_Summary_Report/figures/*", emit: figures, optional: true

    script:
    """
    mkdir -p Master_Summary_Report
    mkdir -p Master_Summary_Report/data
    mkdir -p Master_Summary_Report/tables
    mkdir -p Master_Summary_Report/figures

    quarto render ${qmd} \\
      -P project_name="${project_name}" \\
      -P seurat_rds="${seurat_rds}" \\
      -P export_cells_file="${export_cells}" \\
      -P vdj_qc_per_sample_compact_file="${vdj_qc_per_sample_compact}" \\
      -P vdj_qc_before_after_summary_file="${vdj_qc_before_after_summary}" \\
      -P vdj_qc_sample_sheet_resolved_file="${vdj_qc_sample_sheet_resolved}" \\
      -P vdj_qc_clone_rank_abundance_file="${vdj_qc_clone_rank_abundance}" \\
      -P vdj_qc_before_after_retention_fig="${vdj_qc_before_after_retention_fig}" \\
      -P vdj_qc_pairing_bar_fig="${vdj_qc_pairing_bar_fig}" \\
      -P vdj_qc_clone_rank_abundance_fig="${vdj_qc_clone_rank_abundance_fig}" \\
      -P vdj_qc_multiple_chains_fig="${vdj_qc_multiple_chains_fig}" \\
      -P label_col="${params.label_col}" \\
      -P sample_col="${params.sample_col}" \\
      -P patient_col="${params.patient_col}" \\
      -P condition_col="${params.condition_col}" \\
      -P timepoint_col="${params.timepoint_col}" \\
      -P outdir="Master_Summary_Report"
    """
}