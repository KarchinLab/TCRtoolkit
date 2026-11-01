/*
 * Bridge: REPORT_BUNDLE
 *
 * Gathers the module tables the report notebooks read into one report_bundle/<module>/
 * directory, on every route: the run-level sections (QC & Run Info overview and VDJ QC,
 * repertoire structure, overlap, sequence clustering) and the single-cell sections. Several
 * modules write same-named tables (method_presence_summary.tsv, annotation_*.tsv), which
 * would collide if staged flat into RENDER_NOTEBOOK; a single directory input also keeps
 * that process's inputs fixed regardless of which modules ran. A module that did not run
 * contributes an empty directory, which the notebooks treat as "absent". Only tables are
 * copied (no .rds objects).
 *
 * run_info.tsv records the route ("bulk", "vdj_only" or "full") and the run parameters the
 * notebooks quote.
 */
process REPORT_BUNDLE {
    tag "${params.project_name}"
    label 'process_single'
    container "${params.container}"

    input:
    path cells,             stageAs: 'in/cells/*'
    path vdj_qc_tables,     stageAs: 'in/vdj_qc/*'
    path pseudobulk_tables, stageAs: 'in/pseudobulk/*'
    path sample_tables,     stageAs: 'in/sample/*'
    path compare_tables,    stageAs: 'in/compare/*'
    path repertoire_tables, stageAs: 'in/repertoire/*'
    path tcell_tables,      stageAs: 'in/tcell/*'
    path mergedvdj_tables,  stageAs: 'in/mergedvdj/*'
    path consensus_tables,  stageAs: 'in/consensus/*'
    path rollup_tables,     stageAs: 'in/rollup/*'
    path conga_tables,      stageAs: 'in/conga/*'
    path tcri_tables,       stageAs: 'in/tcri/*'
    val  route

    output:
    path "report_bundle", emit: bundle

    script:
    """
    mkdir -p report_bundle
    for d in cells vdj_qc pseudobulk sample compare repertoire tcell mergedvdj consensus rollup conga tcri; do
        mkdir -p report_bundle/\$d
        for f in in/\$d/*.tsv in/\$d/*.csv in/\$d/*.txt; do
            if [ -s "\$f" ]; then cp -L "\$f" report_bundle/\$d/; fi
        done
    done
    echo "${route}" > report_bundle/route.txt
    printf 'key\tvalue\nroute\t%s\nproject_name\t%s\ntcrdist_radius\t%s\n' "${route}" "${params.project_name}" "${params.tcrdist_radius}" > report_bundle/run_info.tsv
    """
}
