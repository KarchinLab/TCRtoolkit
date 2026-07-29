/*
 * VALIDATE_PARAMS
 *
 * Lightweight required-input checks. The nf-schema plugin's validateParameters() is not
 * used here: the version available in this environment (2.4.2) cannot parse the
 * draft-2020-12 / $defs schema, and its session observer aborts every run. The pipeline's
 * own guards cover the essentials (single-cell input checks live in SINGLECELL_WORKFLOW).
 */
workflow VALIDATE_PARAMS {

    main:
    def mode = (params.mode ?: 'bulk').toLowerCase()

    if (mode == 'bulk' && !params.samplesheet) {
        error "Bulk mode requires --samplesheet <file>. For single-cell input use --mode singlecell."
    }
    if (!params.outdir) {
        error "Please provide an output directory with --outdir."
    }
}
