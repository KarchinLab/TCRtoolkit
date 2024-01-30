
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CALC_COMPARE  } from '../../modules/local/calc_compare.nf'
include { PLOT_COMPARE  } from '../../modules/local/plot_compare.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow COMPARE {

    println("Welcome to the BULK TCRSEQ pipeline! -- COMPARE ")

    take:
    sample_utf8
    meta_data

    main:
    CALC_COMPARE( sample_utf8,
                  meta_data )

    PLOT_COMPARE( CALC_COMPARE.out.jaccard_mat,
                  CALC_COMPARE.out.sorensen_mat,
                  CALC_COMPARE.out.morisita_mat )
    
    // emit:
    // compare_stats_html
    // versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}