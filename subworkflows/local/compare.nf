
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TCRSHARING_CALC } from '../../modules/local/compare/tcrsharing'
include { OLGA_MERGE as TCRSHARING_OLGA_MERGE } from '../../modules/local/olga'
include { CLONAL_DYNAMICS } from '../../modules/local/compare/clonal_dynamics'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow COMPARE {
    take:
    concat_cdr3_sorted
    cdr3_pgen

    main:
    TCRSHARING_OLGA_MERGE (concat_cdr3_sorted, cdr3_pgen)

    TCRSHARING_CALC(
        TCRSHARING_OLGA_MERGE.out.concat_cdr3_pgen
    )

    CLONAL_DYNAMICS(
        concat_cdr3_sorted,
        file(params.samplesheet)
    )

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
}