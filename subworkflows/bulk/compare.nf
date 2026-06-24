include { TCRSHARING_CALC; TCRSHARING_HISTOGRAM; TCRSHARING_SCATTERPLOT } from '../../modules/bulk/compare/tcrsharing'
include { OLGA_MERGE as TCRSHARING_OLGA_MERGE }                          from '../../modules/bulk/olga/main'

workflow COMPARE {
    take:
    concat_cdr3_sorted
    cdr3_pgen

    main:
    TCRSHARING_OLGA_MERGE( concat_cdr3_sorted, cdr3_pgen )

    TCRSHARING_CALC( TCRSHARING_OLGA_MERGE.out.concat_cdr3_pgen )

    TCRSHARING_HISTOGRAM( TCRSHARING_CALC.out.shared_cdr3 )

    TCRSHARING_SCATTERPLOT( TCRSHARING_CALC.out.shared_cdr3 )
}
