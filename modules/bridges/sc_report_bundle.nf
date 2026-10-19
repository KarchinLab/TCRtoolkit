/*
 * Bridge: SC_REPORT_BUNDLE
 *
 * Gathers the single-cell module tables the report notebooks' SC sections read into one
 * sc_bundle/<module>/ directory. Several modules write same-named tables
 * (method_presence_summary.tsv, annotation_*.tsv), which would collide if staged flat
 * into RENDER_NOTEBOOK; a single directory input also keeps that process's inputs
 * fixed regardless of which modules ran. A module that did not run contributes an
 * empty directory, which the notebooks treat as "absent".
 */
process SC_REPORT_BUNDLE {
    tag "${params.project_name}"
    label 'process_single'
    container "${params.container}"

    input:
    path cells,             stageAs: 'in/cells/*'
    path vdj_qc_tables,     stageAs: 'in/vdj_qc/*'
    path pseudobulk_tables, stageAs: 'in/pseudobulk/*'
    path tcell_tables,      stageAs: 'in/tcell/*'
    path mergedvdj_tables,  stageAs: 'in/mergedvdj/*'
    path consensus_tables,  stageAs: 'in/consensus/*'
    path rollup_tables,     stageAs: 'in/rollup/*'
    path conga_tables,      stageAs: 'in/conga/*'
    path tcri_tables,       stageAs: 'in/tcri/*'
    val  route

    output:
    path "sc_bundle", emit: bundle

    script:
    """
    mkdir -p sc_bundle
    for d in cells vdj_qc pseudobulk tcell mergedvdj consensus rollup conga tcri; do
        mkdir -p sc_bundle/\$d
        if [ -d in/\$d ]; then
            find -L in/\$d -maxdepth 1 -type f -size +0 -exec cp -L {} sc_bundle/\$d/ \\;
        fi
    done
    echo "${route}" > sc_bundle/route.txt
    """
}
