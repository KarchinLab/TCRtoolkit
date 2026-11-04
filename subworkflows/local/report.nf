/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { RENDER_NOTEBOOK } from '../../modules/local/report/render_notebook'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Each main notebook ends with {{< include ./_sc_<name>.qmd >}}: the single-cell-only
// section appended to it (CoNGA, TCRi, phenotype, cell-level clustering).
def scSection(notebookName) {
    [
        template_qc              : '_sc_qc.qmd',
        template_details_sample  : '_sc_details_sample.qmd',
        template_details_compare : '_sc_details_compare.qmd',
        template_details_patient : '_sc_details_patient.qmd',
    ][notebookName]
}

// Bulk stand-in for a notebook's single-cell section, staged under the section's name.
def scBulkSection(notebookName) {
    [ template_details_compare : '_bulk_details_compare.qmd' ][notebookName]
}

// Run-level sections (moved from the Master Summary) each notebook includes, all reading
// report_bundle/ through _report_bundle.qmd. template_sample.qmd is the body of the
// sample details report.
def runSections(notebookName) {
    [
        template_qc              : ['_qc_overview.qmd', '_qc_technical.qmd'],
        template_details_sample  : ['_sample_gene_usage.qmd', '_sample_clonal_structure.qmd', '_violin_heatmap.qmd'],
        template_details_compare : ['_compare_overlap.qmd'],
        template_details_patient : ['_patient_clustering.qmd', '_violin_heatmap.qmd'],
    ][notebookName] ?: []
}

// Single-channel listener: receives tuples of
// (notebook_template, files_to_stage, staged_layout)
// and renders each notebook to HTML.
workflow REPORT {

    take:
    ch_reports  // channel of tuples: tuple(path notebook_template, path files_to_stage, staged_layout)
    samplesheet // value channel: samplesheet the notebooks read as sample_table
    bundle      // value channel: REPORT_BUNDLE's report_bundle/ directory (every route)
    route       // "bulk", "vdj_only" or "full"

    main:
    def notebooks_dir = file(params.template_qc).parent
    // Shared look: _report_style.qmd includes the CSS and plot theme partials.
    def style_files = ['_report_css.qmd', '_report_theme.qmd'].collect { f -> file("${notebooks_dir}/${f}") }

    // Quarto includes are static, so the SC section is always included: with SC data the
    // real partial is staged alongside the bundle tables it reads; on the bulk route an
    // empty stub is staged under the partial's name and the notebook renders as before.
    // The Discovery Brief is not published on its own: it is flattened to markdown and
    // becomes the opening sections of the Master Summary.
    ch_with_sc = ch_reports
        .combine(bundle)
        .map { notebook, files, layout, bundle_dir ->
            def name = notebook.baseName
            def extra = style_files + [bundle_dir]
            def run = runSections(name)
            if (run) {
                extra += ([file("${notebooks_dir}/_report_bundle.qmd")] +
                          run.collect { f -> file("${notebooks_dir}/${f}") })
            }
            def section = scSection(name)
            def to_markdown = name == 'template_discovery_brief'
            if (!section) return [notebook, files + extra, layout, to_markdown]
            if (route != 'bulk') {
                def sc_files = [file("${notebooks_dir}/${section}"),
                                file("${notebooks_dir}/_sc_common.qmd")]
                [notebook, files + extra + sc_files, layout, to_markdown]
            } else {
                // Bulk data: the stub, except for the compare report, whose clone-level
                // longitudinal section has a bulk counterpart.
                def stand_in = scBulkSection(name) ?: '_sc_stub.qmd'
                [notebook, files + extra + [file("${notebooks_dir}/${stand_in}")],
                 layout + [[section, stand_in]], to_markdown]
            }
        }

    RENDER_NOTEBOOK(
        ch_with_sc,
        params.project_name,
        workflow.commandLine,
        samplesheet
    )

    emit:
        report_html  = RENDER_NOTEBOOK.out.report_html
        figures      = RENDER_NOTEBOOK.out.figures
        // The Discovery Brief as markdown for the Master Summary (empty without it).
        discovery_md = RENDER_NOTEBOOK.out.markdown
}
