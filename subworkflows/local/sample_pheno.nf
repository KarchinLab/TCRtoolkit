
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMPLE_CALC } from '../../modules/local/sample/sample_calc'
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

    SAMPLE_CALC( sample_map )

    SAMPLE_CALC.out.sample_csv
        .collectFile(name: 'sample_stats.csv', sort: true, 
                     storeDir: "${params.outdir}/sample_phenotype")
        .set { sample_stats_csv }

    SAMPLE_CALC.out.v_family_csv
        .collectFile(name: 'v_family.csv', sort: true,
                     storeDir: "${params.outdir}/sample_phenotype")
        .set { v_family_csv }

    SAMPLE_CALC.out.d_family_csv
        .collectFile(name: 'd_family.csv', sort: true,
                     storeDir: "${params.outdir}/sample_phenotype")
        .set { d_family_csv }

    SAMPLE_CALC.out.j_family_csv
        .collectFile(name: 'j_family.csv', sort: true,
                     storeDir: "${params.outdir}/sample_phenotype")
        .set { j_family_csv }

    /////// =================== SAMPLE NOTEBOOK ===================  ///////

    // SAMPLE_PLOT (
    //     pheno_samplesheet,
    //     file(params.sample_stats_template),
    //     sample_stats_csv,
    //     v_family_csv
    //     )

}