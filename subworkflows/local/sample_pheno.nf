
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMPLE_CALC as SAMPLE_CALC_PHENO } from '../../modules/local/sample/sample_calc'
include { SAMPLE_AGGREGATE as SAMPLE_AGG_STAT_PHENO } from '../../modules/local/sample/sample_aggregate' 
include { SAMPLE_AGGREGATE as SAMPLE_AGG_V_PHENO } from '../../modules/local/sample/sample_aggregate'
include { SAMPLE_AGGREGATE as SAMPLE_AGG_D_PHENO } from '../../modules/local/sample/sample_aggregate'
include { SAMPLE_AGGREGATE as SAMPLE_AGG_J_PHENO } from '../../modules/local/sample/sample_aggregate'
// include { SAMPLE_PLOT } from '../../modules/local/sample/sample_plot'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SAMPLE_PHENO {

    take:
    sample_map
    pheno_samplesheet

    main:

    /////// =================== CALC SAMPLE ===================  ///////

    SAMPLE_CALC_PHENO( sample_map )

    SAMPLE_CALC_PHENO.out.sample_csv.collect().set { sample_csv_files }
    SAMPLE_CALC_PHENO.out.v_family_csv.collect().set { v_family_csv_files }
    SAMPLE_CALC_PHENO.out.d_family_csv.collect().set { d_family_csv_files }
    SAMPLE_CALC_PHENO.out.j_family_csv.collect().set { j_family_csv_files }

    SAMPLE_AGG_STAT_PHENO(sample_csv_files, "sample_stats.csv")
    SAMPLE_AGG_V_PHENO(v_family_csv_files, "v_family.csv")
    SAMPLE_AGG_D_PHENO(d_family_csv_files, "d_family.csv")
    SAMPLE_AGG_J_PHENO(j_family_csv_files, "j_family.csv")


    /////// =================== SAMPLE NOTEBOOK ===================  ///////

    // SAMPLE_PLOT (
    //     pheno_samplesheet,
    //     file(params.sample_stats_template),
    //     sample_stats_csv,
    //     v_family_csv
    //     )

}