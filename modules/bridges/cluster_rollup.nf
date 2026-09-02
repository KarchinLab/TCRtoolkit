/*
 * CLUSTER_ROLLUP
 *
 * GIANA, GLIPH2 and TCRdist3 write raw per-patient / per-sample output but no rollup
 * table, so MASTER_SUMMARY's giana_summary_file / gliph2_summary_file /
 * tcrdist3_summary_file had nothing to read and those modules always reported as
 * absent. This process summarises them into the filenames the .qmd expects.
 *
 * Every output is `optional: true`: a route where a method did not run simply yields
 * no file, and the workflow substitutes the NO_FILE placeholder that the .qmd already
 * treats as missing.
 */
process CLUSTER_ROLLUP {
    tag "${project_name}"
    label 'process_single'
    container "${params.container}"

    publishDir "${params.outdir}/bridge/cluster_rollup", mode: 'copy', overwrite: true

    input:
      path giana_files,   stageAs: 'giana/*'
      path gliph2_files,  stageAs: 'gliph2/*'
      path tcrdist_files, stageAs: 'tcrdist/*'
      path tcrdist_mats,  stageAs: 'tcrdistmat/*'
      path export_cells,  stageAs: 'in_export_cells.tsv'
      // Staged rather than invoked as ${projectDir}/bin/... so its contents take part in
      // the task hash: otherwise editing the script leaves -resume silently reusing the
      // previous output.
      path rollup_script, stageAs: 'cluster_rollup.py'
      val  project_name

    output:
      path "giana_summary_rollup.tsv",        emit: giana_summary,        optional: true
      path "gliph2_summary_rollup.tsv",       emit: gliph2_summary,       optional: true
      path "tcrdist3_summary_rollup.tsv",     emit: tcrdist3_summary,     optional: true
      path "method_presence_summary.tsv",     emit: method_presence,      optional: true
      path "method_cluster_counts.tsv",       emit: method_cluster_counts, optional: true
      path "annotation_giana_summary.tsv",    emit: annotation_giana,     optional: true
      path "annotation_gliph2_summary.tsv",   emit: annotation_gliph2,    optional: true
      path "annotation_tcrdist3_summary.tsv", emit: annotation_tcrdist3,  optional: true
      path "giana_cluster_detail.tsv",        emit: giana_detail,        optional: true
      path "gliph2_motif_detail.tsv",         emit: gliph2_detail,       optional: true
      path "tcrdist_cluster_detail.tsv",      emit: tcrdist_detail,      optional: true
      path "giana_vgene_usage.tsv",           emit: giana_vgene,         optional: true
      path "gliph2_vgene_usage.tsv",          emit: gliph2_vgene,        optional: true

    script:
    // nullglob so an empty staging directory expands to nothing rather than a literal
    // "giana/*". Guards are written as if-blocks, not `test && assign`: Nextflow runs the
    // script under `set -e`, where a failing `&&` chain would abort the whole task.
    """
    shopt -s nullglob

    giana_args=()
    gliph_args=()
    tcrdist_args=()
    matrix_args=()

    g=( giana/* )
    if [ \${#g[@]} -gt 0 ]; then giana_args=( --giana "\${g[@]}" ); fi

    p=( gliph2/* )
    if [ \${#p[@]} -gt 0 ]; then gliph_args=( --gliph2 "\${p[@]}" ); fi

    t=( tcrdist/* )
    if [ \${#t[@]} -gt 0 ]; then tcrdist_args=( --tcrdist "\${t[@]}" ); fi

    m=( tcrdistmat/* )
    if [ \${#m[@]} -gt 0 ]; then matrix_args=( --tcrdist-matrix "\${m[@]}" ); fi

    python3 cluster_rollup.py \\
        --outdir . \\
        --export-cells in_export_cells.tsv \\
        --tcrdist-radius ${params.tcrdist_radius} \\
        \${matrix_args[@]+"\${matrix_args[@]}"} \\
        \${giana_args[@]+"\${giana_args[@]}"} \\
        \${gliph_args[@]+"\${gliph_args[@]}"} \\
        \${tcrdist_args[@]+"\${tcrdist_args[@]}"}
    """
}
