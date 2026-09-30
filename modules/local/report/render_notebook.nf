// Generic process to render a Quarto notebook to HTML
process RENDER_NOTEBOOK {
    tag "${notebook.getBaseName()}"

    // template_details_compare and template_discovery_brief process full VDJdb match
    // files (up to ~1GB/sample) and build per-patient network graphs; the rest
    // aggregate small summary tables and finish in under 5 minutes.
    cpus   { ['template_details_compare', 'template_discovery_brief'].contains(notebook.baseName) ? 16 * task.attempt : 4 * task.attempt }
    memory { ['template_details_compare', 'template_discovery_brief'].contains(notebook.baseName) ? 256.GB * task.attempt : 16.GB * task.attempt }
    time   { ['template_details_compare', 'template_discovery_brief'].contains(notebook.baseName) ? 16.h : 4.h }

    input:
    // path(files) stages files flat in the root dir; staged_layout optionally
    // symlinks them into a project_dir-style subdirectory tree for notebooks that
    // read from a nested project_dir/subdir/file layout instead of bare filenames.
    // It's a list of [dest_path, source_basename] pairs (e.g.
    // ["sample/sample_stats.csv", "sample_stats.csv"]) rather than just a dest
    // path, because source and dest basenames can legitimately differ - e.g.
    // template_discovery_brief.qmd includes a generic template_pheno.qmd, which
    // is resolved to the real notebook that applies (template_pheno_bulk.qmd)
    // and symlinked to that shared destination name.
    tuple path(notebook), path(files), val(staged_layout)
    val project_name
    val workflow_cmd
    // Fixed name avoids colliding with a same-named file in `files`.
    path samplesheet, stageAs: 'render_notebook_samplesheet.csv'

    output:
    path "${notebook.getBaseName()}.html", emit: report_html

    script:
    // Joined with '\n    ' (not just '\n') so each generated line lands at the
    // same 4-space indent as the surrounding script block below, keeping the
    // rendered task script readable in .command.sh.
    def stage_cmds = staged_layout.collect { dest, src ->
        "mkdir -p \"\$(dirname '${dest}')\"; ln -sf \"\$PWD/${src}\" '${dest}'"
    }.join('\n    ')
    def project_dir_arg = staged_layout ? "-P project_dir:'.'" : ''
    """
    ${stage_cmds}
    ## render qmd report to html
    quarto render ${notebook} \\
        -P project_name:${project_name} \\
        -P workflow_cmd:'${workflow_cmd}' \\
        -P sample_table:${samplesheet} \\
        -P subject_col:'${params.subject_col}' \\
        -P timepoint_col:'${params.timepoint_col}' \\
        -P timepoint_order_col:'${params.timepoint_order_col}' \\
        -P timepoint_order:'${params.timepoint_order}' \\
        -P alias_col:'${params.alias_col}' \\
        ${project_dir_arg} \\
        --to html
    """

    stub:
    """
    touch ${notebook.getBaseName()}.html
    """
}
