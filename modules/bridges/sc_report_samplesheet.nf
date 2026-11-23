/*
 * Bridge: SC_REPORT_SAMPLESHEET
 *
 * Builds the bulk-style samplesheet (sample, alias, patient, timepoint, origin, file) the
 * report notebooks read, one row per QC-passing pseudobulk unit, with metadata taken from
 * the user's SC sample sheet.
 */
process SC_REPORT_SAMPLESHEET {
    tag "report samplesheet"
    label 'process_single'
    container "${params.container}"
    publishDir "${params.outdir}/bridge/report_samplesheet", mode: 'copy', overwrite: true

    input:
    path units
    path sample_sheet, stageAs: 'sc_sample_sheet.csv'
    path cells,        stageAs: 'in_export_cells.tsv'

    output:
    path "report_samplesheet.csv", emit: samplesheet

    script:
    def ss = (sample_sheet.size() > 0) ? "--sample-sheet ${sample_sheet}" : ""
    def ce = (cells.size() > 0) ? "--cells ${cells}" : ""
    """
    sc_report_samplesheet.py ${units} ${ss} ${ce} --timepoint-col '${params.timepoint_col}' --out report_samplesheet.csv
    """
}
