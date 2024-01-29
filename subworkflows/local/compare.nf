
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
    

    
    /////// =================== PLOT COMPARE ===================  ///////
    // PLOT_COMPARE(
    //     file(params.sample_table),
    //     sample_stats_csv
    //     )
    
    // emit:
    // sample_stats_csv
    // v_family_csv
    // sample_meta_csv
    // versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}