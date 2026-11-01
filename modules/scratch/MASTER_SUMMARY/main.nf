process MASTER_SUMMARY {
    tag "${project_name}"
    label 'process_medium'
    container "${params.container}"

    // Published with the report notebooks (report/<name>/), as report/master_summary/.
    publishDir "${params.outdir}/report/master_summary", mode: 'copy', overwrite: true

    input:
      // Core per-cell inputs. stageAs with a real extension matters: the .qmd picks its
      // reader from the extension, and knitr/quarto infer image type the same way.
      path seurat_rds,        stageAs: 'in_seurat_rds.rds'
      path export_cells,      stageAs: 'in_export_cells.tsv'

      // Each upstream module contributes its whole tables/ directory under
      // intables/<module>/: a module that did not run on this route simply stages nothing,
      // the .qmd finds no files, and the affected sections render an explanatory note. Only
      // the modules behind the sections still rendered here (consensus clustering and the
      // GEX section); the run-level QC, repertoire and clustering tables go to the report
      // notebooks through REPORT_BUNDLE.
      path rollup_tables,     stageAs: 'intables/rollup/*'
      path tcell_tables,      stageAs: 'intables/tcell/*'
      path conga_tables,      stageAs: 'intables/conga/*'
      path consensus_tables,  stageAs: 'intables/consensus/*'
      path tcri_tables,       stageAs: 'intables/tcri/*'
      path clonal_tables,     stageAs: 'intables/clonal/*'

      // The Discovery Brief flattened to markdown by RENDER_NOTEBOOK
      // (bin/notebook_to_markdown.py), or NO_FILE without the report notebooks.
      path discovery_md,      stageAs: 'in_discovery_brief.md'

      // Master_Summary_Report.qmd plus its partials (_ms_setup, _sc_gex, _sc_gex_stub) and
      // the standalone SC_GEX_Report.qmd wrapper.
      path qmd_files,         stageAs: 'src/*'
      val  barrier_done
      val  project_name

    output:
      path "Master_Summary_Report.html", emit: report_html
      path "Master_Summary_Report/tables/*",  emit: tables,  optional: true
      path "Master_Summary_Report/figures/*", emit: figures, optional: true
      // GEX route only: the single-cell GEX section rendered on its own.
      path "SC_GEX_Report.html",              emit: gex_html, optional: true
      path "SC_GEX_Report/tables/*",          emit: gex_tables, optional: true

    script:
    """
    mkdir -p Master_Summary_Report/tables Master_Summary_Report/figures
    cp -L src/*.qmd .

    # The GEX section is a static include: the T-cell integration export only exists when a
    # GEX object was supplied, so without it the stub stands in under the partial's name.
    gex=false
    if [ -s intables/tcell/tcr_export_cells_with_embedding.tsv ]; then gex=true; else cp _sc_gex_stub.qmd _sc_gex.qmd; fi

    # The Discovery Brief's sections open the report; without the report notebooks a note
    # stands in under the included name.
    if [ -s in_discovery_brief.md ]; then
      cp -L in_discovery_brief.md _discovery_brief.md
    else
      printf '%s\\n' '::: {.callout-note appearance="simple"}' \
        'The repertoire findings (clonality, overlap, clonal dynamics, TCRdist3 networks, VDJdb and publicity) come from the report notebooks, which did not run (`--run_reports false`).' \
        ':::' > _discovery_brief.md
    fi

    common=(
      -P project_name="${project_name}"
      -P modality="${(params.mode ?: 'bulk').toLowerCase()}"
      -P tables_dir="intables"
      -P seurat_rds="${seurat_rds}"
      -P export_cells_file="${export_cells}"
      -P label_col="${params.label_col}"
      -P sample_col="${params.sample_col}"
      -P patient_col="${params.patient_col}"
      -P condition_col="${params.condition_col}"
      -P timepoint_col="${params.timepoint_col}"
      -P timepoint_order="${params.timepoint_order ?: ''}"
      -P tcrdist_radius=${params.tcrdist_radius}
    )

    quarto render Master_Summary_Report.qmd "\${common[@]}" -P outdir="Master_Summary_Report"

    if [ "\$gex" = true ]; then
      mkdir -p SC_GEX_Report/tables SC_GEX_Report/figures
      quarto render SC_GEX_Report.qmd "\${common[@]}" -P outdir="SC_GEX_Report"
    fi
    """
}
