// Generic process to render a Quarto notebook to HTML
process RENDER_NOTEBOOK {
    tag "${notebook.getBaseName()}"

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
    // to_markdown: render, then flatten the executed notebook to <name>.md for another
    // report to include, instead of publishing this one's HTML (the Discovery Brief, shown
    // inside the Master Summary).
    tuple path(notebook), path(files), val(staged_layout), val(to_markdown)
    val project_name
    val workflow_cmd
    // Fixed name avoids colliding with a same-named file in `files`.
    path samplesheet, stageAs: 'render_notebook_samplesheet.csv'

    output:
    path "${notebook.getBaseName()}.html", emit: report_html, optional: true
    path "${notebook.getBaseName()}.md",   emit: markdown,    optional: true
    // Notebooks that save figures write them under ./figures/ (PNGs for static
    // plots, standalone interactive .html fragments for Plotly figures kept
    // out of the main report document). Optional: not every notebook does.
    path "figures", emit: figures, optional: true

    script:
    // Joined with '\n    ' (not just '\n') so each generated line lands at the
    // same 4-space indent as the surrounding script block below, keeping the
    // rendered task script readable in .command.sh.
    def stage_cmds = staged_layout.collect { dest, src ->
        "mkdir -p \"\$(dirname '${dest}')\"; ln -sf \"\$PWD/${src}\" '${dest}'"
    }.join('\n    ')
    def project_dir_arg = staged_layout ? "-P project_dir:'.'" : ''
    def base = notebook.getBaseName()
    def keep_ipynb_arg = to_markdown ? '-M keep-ipynb:true' : ''
    def flatten_cmd = to_markdown
        ? "notebook_to_markdown.py ${base}.quarto_ipynb ${base}.md && rm -f ${base}.html"
        : ''
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
        ${keep_ipynb_arg} \\
        --to html
    ${flatten_cmd}
    """

    stub:
    """
    touch ${notebook.getBaseName()}.html
    """
}
