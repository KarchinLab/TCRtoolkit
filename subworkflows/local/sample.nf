
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMPLE_CALC; SAMPLE_CALC_PIVOT } from '../../modules/local/sample/sample_calc'
include { OLGA_SAMPLE_MERGE } from '../../modules/local/olga'
include { CONVERGENCE } from '../../modules/local/sample/convergence'
include { TCRPHENO } from '../../modules/local/sample/tcrpheno'
include { VDJDB_GET; VDJDB_VDJMATCH } from '../../modules/local/sample/tcrspecificity'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SAMPLE {

    take:
    processed_samples
    per_sample_stats
    cdr3_pgen
    olga_stats
    run_tcrpheno    // bool: false on the single-cell VDJ+GEX route, which has measured phenotypes

    main:

    /////// =================== CALC SAMPLE ===================  ///////

    def sample_calc_input = processed_samples.join(per_sample_stats)
    SAMPLE_CALC( sample_calc_input )

    // sample_stats.csv is published by REPERTOIRE_METRICS (bulktcr_analysis.nf), which
    // adds the diversity columns to these per-sample rows.

    def v_family_agg = SAMPLE_CALC.out.v_family_csv
        .collectFile(name: "v_family_long.csv", keepHeader: true, skip: 1, sort: true)

    def d_family_agg = SAMPLE_CALC.out.d_family_csv
        .collectFile(name: "d_family_long.csv", keepHeader: true, skip: 1, sort: true)

    def j_family_agg = SAMPLE_CALC.out.j_family_csv
        .collectFile(name: "j_family_long.csv", keepHeader: true, skip: 1, sort: true)

    SAMPLE_CALC_PIVOT( v_family_agg, d_family_agg, j_family_agg )

    // TCRdist3 runs per patient (PATIENT subworkflow), which also writes each sample's
    // block of the patient matrix for the per-sample report sections.

    OLGA_SAMPLE_MERGE ( processed_samples,
        cdr3_pgen.first(),
        olga_stats )

    CONVERGENCE ( processed_samples )

    if (run_tcrpheno) {
        TCRPHENO ( processed_samples )
    }

    VDJDB_GET ()

    VDJDB_VDJMATCH (processed_samples, VDJDB_GET.out.ref_db)

    emit:
        sample_csv        = SAMPLE_CALC.out.sample_csv
        v_family          = SAMPLE_CALC_PIVOT.out.v_family_wide
        j_family          = SAMPLE_CALC_PIVOT.out.j_family_wide
        olga_files        = OLGA_SAMPLE_MERGE.out.olga_pgen.map { _sample_meta, pgen -> pgen }
                                .collect()
        vdjdb_files       = VDJDB_VDJMATCH.out.vdjmatch_txt
                                .mix( VDJDB_VDJMATCH.out.annot_summary )
                                .collect()
        convergence_files = CONVERGENCE.out.convergence_output.collect()
        tcrpheno_files    = run_tcrpheno ? TCRPHENO.out.tcrpheno_output.map { _sample_meta, f -> f }.collect() : channel.value([])
}