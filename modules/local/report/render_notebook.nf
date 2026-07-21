// Generic process to render a Quarto notebook to HTML
process RENDER_NOTEBOOK {
    tag "${notebook.getBaseName()}"
    label 'process_single'

    input:
    // path(files) stages files flat in the root dir; staged_layout optionally
    // symlinks them into a project_dir-style subdirectory tree for notebooks that
    // read from a nested project_dir/subdir/file layout instead of bare filenames.
    // It's a list of [dest_path, source_basename] pairs (e.g.
    // ["sample/sample_stats.csv", "sample_stats.csv"]) rather than just a dest
    // path, because source and dest basenames can legitimately differ - e.g.
    // per-patient files staged under a shared basename (like GLIPH2's
    // "patientA_all_motifs.txt") get symlinked to their unprefixed, per-patient
    // subdirectory destination ("gliph2/patientA/all_motifs.txt").
    tuple path(notebook), path(files), val(staged_layout)
    val project_name
    val workflow_cmd

    output:
    path "${notebook.getBaseName()}.html", emit: report_html

    script:
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
        -P sample_table:${file(params.samplesheet)} \\
        -P subject_col:'${params.subject_col}' \\
        -P timepoint_col:'${params.timepoint_col}' \\
        -P timepoint_order_col:'${params.timepoint_order_col}' \\
        -P alias_col:'${params.alias_col}' \\
        ${project_dir_arg} \\
        --to html
    """

    stub:
    """
    touch ${notebook.getBaseName()}.html
    """
}
