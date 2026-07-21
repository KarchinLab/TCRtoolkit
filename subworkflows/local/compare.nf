
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { TCRSHARING_CALC; TCRSHARING_HISTOGRAM; TCRSHARING_SCATTERPLOT } from '../../modules/local/compare/tcrsharing'
include { OLGA_MERGE as TCRSHARING_OLGA_MERGE } from '../../modules/local/olga'

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

    TCRSHARING_HISTOGRAM(
        TCRSHARING_CALC.out.shared_cdr3
    )

    TCRSHARING_SCATTERPLOT(
        TCRSHARING_CALC.out.shared_cdr3
    )

    emit:
        shared_cdr3 = TCRSHARING_CALC.out.shared_cdr3
}