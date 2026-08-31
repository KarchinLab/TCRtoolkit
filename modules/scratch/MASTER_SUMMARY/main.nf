process MASTER_SUMMARY {
    tag "${project_name}"
    label 'process_medium'
    container "${params.container}"

    publishDir "${params.outdir}/Master_Summary", mode: 'copy', overwrite: true

    input:
      // Core per-cell inputs. stageAs with a real extension matters: the .qmd picks its
      // reader from the extension, and knitr/quarto infer image type the same way.
      path seurat_rds,        stageAs: 'in_seurat_rds.rds'
      path export_cells,      stageAs: 'in_export_cells.tsv'

      // Each upstream module contributes its whole tables/ directory under
      // intables/<module>/. This replaces the previous ~30 individually-wired file params:
      // a module that did not run on this route simply stages nothing, the .qmd finds no
      // files, and the affected sections render an explanatory note. Adding a new upstream
      // table no longer requires another process input.
      path vdj_qc_tables,     stageAs: 'intables/vdj_qc/*'
      path pseudobulk_tables, stageAs: 'intables/pseudobulk/*'
      path sample_tables,     stageAs: 'intables/sample/*'
      path compare_tables,    stageAs: 'intables/compare/*'
      path rollup_tables,     stageAs: 'intables/rollup/*'
      path repertoire_tables, stageAs: 'intables/repertoire/*'
      path tcell_tables,      stageAs: 'intables/tcell/*'
      path conga_tables,      stageAs: 'intables/conga/*'
      path consensus_tables,  stageAs: 'intables/consensus/*'
      path tcri_tables,       stageAs: 'intables/tcri/*'

      path qmd
      val  barrier_done
      val  project_name

    output:
      path "Master_Summary_Report.html", emit: report_html
      path "Master_Summary_Report/tables/*",  emit: tables,  optional: true
      path "Master_Summary_Report/figures/*", emit: figures, optional: true

    script:
    """
    mkdir -p Master_Summary_Report/tables Master_Summary_Report/figures

    quarto render ${qmd} \\
      -P project_name="${project_name}" \\
      -P tables_dir="intables" \\
      -P seurat_rds="${seurat_rds}" \\
      -P export_cells_file="${export_cells}" \\
      -P outdir="Master_Summary_Report" \\
      -P label_col="${params.label_col}" \\
      -P sample_col="${params.sample_col}" \\
      -P patient_col="${params.patient_col}" \\
      -P condition_col="${params.condition_col}" \\
      -P timepoint_col="${params.timepoint_col}" \\
      -P tcrdist_radius=${params.tcrdist_radius}
    """
}
