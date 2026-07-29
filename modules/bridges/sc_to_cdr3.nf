/*
 * Bridge 1b: SC_TO_CDR3
 *
 * Converts TCELL_INTEGRATION's per-cell export_cells.tsv directly into
 * per-sample CDR3b/TRBV/TRBJ/counts/sample TSVs and a concatenated file,
 * bypassing the ANNOTATE_PROCESS column-rename step.
 *
 * Replaces the SC_TO_BULK → ANNOTATE_PROCESS chain in full-SC mode.
 */

process SC_TO_CDR3 {
    tag "SC → CDR3 (Bridge 1b)"
    label 'process_low'
    container "${params.container}"
    publishDir "${params.outdir}/bridge/sc_to_cdr3", mode: 'copy', overwrite: true

    input:
    path export_cells

    output:
    path "per_sample/*_cdr3.tsv", emit: per_sample
    path "concat_cdr3.tsv",       emit: concat_cdr3
    path "unit_map.csv",          emit: unit_map

    script:
    def by_pheno = params.pseudobulk_by_phenotype ? "--by-phenotype" : ""
    """
    sc_to_cdr3.py "${export_cells}" ${by_pheno}
    """
}
