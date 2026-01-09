
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COMPARE_CALC  } from '../../modules/local/compare/compare_calc'
include { COMPARE_PLOT  } from '../../modules/local/compare/compare_plot'
include { COMPARE_CONCATENATE  } from '../../modules/local/compare/compare_concatenate'

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
    COMPARE_CALC( samplesheet_resolved,
                    all_sample_files )

    COMPARE_PLOT( samplesheet_resolved,
                  COMPARE_CALC.out.jaccard_mat,
                  COMPARE_CALC.out.sorensen_mat,
                  COMPARE_CALC.out.morisita_mat,
                  file(params.compare_stats_template),
                  params.project_name,
                  all_sample_files
                  )

    COMPARE_CONCATENATE( samplesheet_resolved,
        all_sample_files )

}