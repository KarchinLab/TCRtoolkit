
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TCRSHARING_CALC } from '../../modules/local/compare/tcrsharing'
include { OLGA_MERGE as TCRSHARING_OLGA_MERGE } from '../../modules/local/olga'
include { CLONAL_DYNAMICS } from '../../modules/local/compare/clonal_dynamics'
include { COMPARTMENT_OVERLAP } from '../../modules/local/compare/compartment_overlap'
include { REFERENCE_MATCH } from '../../modules/local/compare/reference_match'

// A reference samplesheet lists its tables in a `file` column, relative to the sheet or
// absolute; anything else is a single reference table and stages on its own.
def isReferenceSheet(sheet) {
    sheet.readLines().first().replace('\uFEFF', '').split(',')*.trim().contains('file')
}

def resolveReference(sheet, String path) {
    def p = path.trim()
    file((p.startsWith('/') || p.contains('://')) ? p : sheet.parent.resolve(p).toString(), checkIfExists: true)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow COMPARE {
    take:
    concat_cdr3_sorted
    cdr3_pgen
    samplesheet     // bulk-style samplesheet (patient / timepoint / origin per sample)

    main:
    TCRSHARING_OLGA_MERGE (concat_cdr3_sorted, cdr3_pgen)

    TCRSHARING_CALC(
        TCRSHARING_OLGA_MERGE.out.concat_cdr3_pgen
    )

    CLONAL_DYNAMICS(
        concat_cdr3_sorted,
        samplesheet
    )

    // Clonotype sharing between each patient's compartments (e.g. tumor / blood). Writes
    // header-only tables and a status line when there is nothing to compare.
    COMPARTMENT_OVERLAP(
        concat_cdr3_sorted,
        samplesheet
    )

    // Optional cross-study comparison against external reference clonotype sets.
    def ch_reference_tables = channel.value([])
    if (params.reference_samplesheet) {
        def ref_sheet = file(params.reference_samplesheet, checkIfExists: true)
        def ref_files = isReferenceSheet(ref_sheet)
            ? channel.fromPath(ref_sheet)
                .splitCsv(header: true)
                .map { row -> resolveReference(ref_sheet, row.file) }
                .collect()
            : channel.value([])
        REFERENCE_MATCH(
            concat_cdr3_sorted,
            samplesheet,
            ref_sheet,
            ref_files
        )
        ch_reference_tables = REFERENCE_MATCH.out.matches
            .mix(REFERENCE_MATCH.out.summary)
            .collect()
    }

    emit:
        shared_cdr3 = TCRSHARING_CALC.out.shared_cdr3
        clonal_dynamics_files = CLONAL_DYNAMICS.out.clonotypes
            .mix(CLONAL_DYNAMICS.out.detailed_expansion)
            .mix(CLONAL_DYNAMICS.out.detailed_contraction)
            .mix(CLONAL_DYNAMICS.out.expanded_summary)
            .mix(CLONAL_DYNAMICS.out.contracting_summary)
            .mix(CLONAL_DYNAMICS.out.expansion_counts)
            .mix(CLONAL_DYNAMICS.out.contraction_counts)
            .collect()
        compartment_tables = COMPARTMENT_OVERLAP.out.summary
            .mix(COMPARTMENT_OVERLAP.out.flow)
            .mix(COMPARTMENT_OVERLAP.out.shared_clones)
            .mix(COMPARTMENT_OVERLAP.out.status)
            .mix(COMPARTMENT_OVERLAP.out.samples)
            .collect()
        reference_tables = ch_reference_tables
}