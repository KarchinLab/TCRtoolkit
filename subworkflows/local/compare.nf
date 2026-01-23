
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { COMPARE_CALC  } from '../../modules/local/compare/compare_calc'
include { COMPARE_PLOT  } from '../../modules/local/compare/compare_plot'
include { COMPARE_CONCATENATE  } from '../../modules/local/compare/compare_concatenate'
include { TCRSHARING_CALC; TCRSHARING_HISTOGRAM; TCRSHARING_SCATTERPLOT } from '../../modules/local/compare/tcrsharing'
include { GLIPH2_TURBOGLIPH; GLIPH2_PLOT } from '../../modules/local/compare/gliph2'
include { GIANA_CALC    } from '../../modules/local/compare/giana'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow COMPARE {
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

    GIANA_CALC(
        COMPARE_CONCATENATE.out.concat_cdr3,
        params.threshold,
        params.threshold_score,
        params.threshold_vgene
    )

    if(params.use_gliph2) {
        GLIPH2_TURBOGLIPH(
            COMPARE_CONCATENATE.out.concat_cdr3
        )
    }

    TCRSHARING_CALC(
        COMPARE_CONCATENATE.out.concat_cdr3
    )

    TCRSHARING_HISTOGRAM(
        TCRSHARING_CALC.out.shared_cdr3
    )

    TCRSHARING_SCATTERPLOT(
        TCRSHARING_CALC.out.shared_cdr3
    )
}