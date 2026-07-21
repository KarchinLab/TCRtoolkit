
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMPLE_CALC; SAMPLE_CALC_PIVOT } from '../../modules/local/sample/sample_calc'
include { SAMPLE_PLOT } from '../../modules/local/sample/sample_plot'
include { TCRDIST3_MATRIX; TCRDIST3_HISTOGRAM_CALC; TCRDIST3_HISTOGRAM_PLOT} from '../../modules/local/sample/tcrdist3'
include { OLGA_SAMPLE_MERGE; OLGA_HISTOGRAM_CALC; OLGA_HISTOGRAM_PLOT } from '../../modules/local/olga'
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

    main:

    /////// =================== CALC SAMPLE ===================  ///////

    def sample_calc_input = processed_samples.join(per_sample_stats)
    SAMPLE_CALC( sample_calc_input )

    def sample_stats_agg = SAMPLE_CALC.out.sample_csv
        .collectFile(name: "sample_stats.csv", keepHeader: true, skip: 1, sort: true, storeDir: "${params.outdir}/sample")

    def v_family_agg = SAMPLE_CALC.out.v_family_csv
        .collectFile(name: "v_family_long.csv", keepHeader: true, skip: 1, sort: true)

    def d_family_agg = SAMPLE_CALC.out.d_family_csv
        .collectFile(name: "d_family_long.csv", keepHeader: true, skip: 1, sort: true)

    def j_family_agg = SAMPLE_CALC.out.j_family_csv
        .collectFile(name: "j_family_long.csv", keepHeader: true, skip: 1, sort: true)

    SAMPLE_CALC_PIVOT( v_family_agg, d_family_agg, j_family_agg )

    /////// =================== PLOT SAMPLE ===================  ///////

    SAMPLE_PLOT (
        file(params.samplesheet),
        file(params.sample_stats_template),
        sample_stats_agg,
        SAMPLE_CALC_PIVOT.out.v_family_wide
        )

    TCRDIST3_MATRIX(
        processed_samples,
        params.matrix_sparsity,
        params.distance_metric,
        file(params.db_path)
    )

    def global_x_max_value = TCRDIST3_MATRIX.out.max_matrix_value
        .map { tcrdist_xmax -> tcrdist_xmax.text.trim().toDouble() }
        .collect()
        .map { values -> values.max() }
    global_x_max_value.view { global_xmax -> "Global x max matrix value: $global_xmax" }

    TCRDIST3_HISTOGRAM_CALC( 
        TCRDIST3_MATRIX.out.tcrdist_output,
        params.matrix_sparsity,
        params.distance_metric,
        global_x_max_value
    )

    def global_y_max_value = TCRDIST3_HISTOGRAM_CALC.out.max_histogram_count
        .map { tcrdist_ymax -> tcrdist_ymax.text.trim().toDouble() }
        .collect()
        .map { values -> values.max() }
    global_y_max_value.view { global_ymax -> "Global y max matrix value: $global_ymax" }

    TCRDIST3_HISTOGRAM_PLOT( 
        TCRDIST3_HISTOGRAM_CALC.out.histogram_data,
        global_y_max_value
    )

    OLGA_SAMPLE_MERGE ( processed_samples,
        cdr3_pgen.first(),
        olga_stats )

    OLGA_HISTOGRAM_CALC ( OLGA_SAMPLE_MERGE.out.olga_pgen, olga_stats )

    def olga_y_max_value = OLGA_HISTOGRAM_CALC.out.olga_ymax
        .map { ymax -> ymax.text.trim().toDouble() }
        .collect()
        .map { values -> values.max() }
    olga_y_max_value.view { olga_ymax -> "Olga y max matrix value: $olga_ymax" }

    OLGA_HISTOGRAM_PLOT( OLGA_HISTOGRAM_CALC.out.olga_histogram, olga_y_max_value )

    CONVERGENCE ( processed_samples )

    TCRPHENO ( processed_samples )

    VDJDB_GET ()

    VDJDB_VDJMATCH (processed_samples, VDJDB_GET.out.ref_db)

    emit:
        sample_csv        = SAMPLE_CALC.out.sample_csv
        v_family          = SAMPLE_CALC_PIVOT.out.v_family_wide
        j_family          = SAMPLE_CALC_PIVOT.out.j_family_wide
        tcrdist_files     = TCRDIST3_MATRIX.out.tcrdist_output.map { _sample_meta, dist -> dist }
                                .mix( TCRDIST3_MATRIX.out.clone_df )
                                .flatten()
                                .collect()
        olga_files        = OLGA_SAMPLE_MERGE.out.olga_pgen.map { _sample_meta, pgen -> pgen }
                                .collect()
        vdjdb_files       = VDJDB_VDJMATCH.out.vdjmatch_txt
                                .mix( VDJDB_VDJMATCH.out.annot_summary )
                                .collect()
        convergence_files = CONVERGENCE.out.convergence_output.collect()
        tcrpheno_files    = TCRPHENO.out.tcrpheno_output.map { _sample_meta, f -> f }.collect()
}