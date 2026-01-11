
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COMPARE_CALC as COMPARE_CALC_PHENO } from '../../modules/local/compare/compare_calc'
// include { COMPARE_PLOT  } from '../../modules/local/compare/compare_plot'
include { COMPARE_CONCATENATE as COMPARE_CONCATENATE_PHENO } from '../../modules/local/compare/compare_concatenate'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow COMPARE_PHENO {

    // println("Welcome to the BULK TCRSEQ pipeline! -- COMPARE ")

    take:
    samplesheet_resolved
    all_sample_files

    main:

    COMPARE_CALC_PHENO( samplesheet_resolved,
                    all_sample_files )

    COMPARE_CONCATENATE_PHENO( samplesheet_resolved,
        all_sample_files )


    /////// =================== COMPARE NOTEBOOK ===================  ///////

    // COMPARE_PLOT( samplesheet_resolved,
    //               COMPARE_CALC_PHENO.out.jaccard_mat,
    //               COMPARE_CALC_PHENO.out.sorensen_mat,
    //               COMPARE_CALC_PHENO.out.morisita_mat,
    //               file(params.compare_stats_template),
    //               params.project_name,
    //               all_sample_files
    //               )

}