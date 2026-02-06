
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SAMPLE_CALC } from '../../modules/local/sample/sample_calc'
include { SAMPLE_PLOT } from '../../modules/local/sample/sample_plot'
include { SAMPLE_AGGREGATE as SAMPLE_AGG_STAT } from '../../modules/local/sample/sample_aggregate' 
include { SAMPLE_AGGREGATE as SAMPLE_AGG_V } from '../../modules/local/sample/sample_aggregate'
include { SAMPLE_AGGREGATE as SAMPLE_AGG_D } from '../../modules/local/sample/sample_aggregate'
include { SAMPLE_AGGREGATE as SAMPLE_AGG_J } from '../../modules/local/sample/sample_aggregate'
include { TCRDIST3_MATRIX; TCRDIST3_HISTOGRAM_CALC; TCRDIST3_HISTOGRAM_PLOT} from '../../modules/local/sample/tcrdist3'
include { OLGA_SAMPLE_MERGE; OLGA_HISTOGRAM_CALC; OLGA_HISTOGRAM_PLOT; OLGA_WRITE_MAX } from '../../modules/local/olga'
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
    sample_map
    cdr3_pgen

    main:

    /////// =================== CALC SAMPLE ===================  ///////

    SAMPLE_CALC( sample_map )

    SAMPLE_CALC.out.sample_csv.collect().set { sample_csv_files }
    SAMPLE_CALC.out.v_family_csv.collect().set { v_family_csv_files }
    SAMPLE_CALC.out.d_family_csv.collect().set { d_family_csv_files }
    SAMPLE_CALC.out.j_family_csv.collect().set { j_family_csv_files }

    SAMPLE_AGG_STAT(sample_csv_files, "sample_stats.csv")
    SAMPLE_AGG_V(v_family_csv_files, "v_family.csv")
    SAMPLE_AGG_D(d_family_csv_files, "d_family.csv")
    SAMPLE_AGG_J(j_family_csv_files, "j_family.csv")

    /////// =================== PLOT SAMPLE ===================  ///////

    SAMPLE_PLOT (
        file(params.samplesheet),
        file(params.sample_stats_template),
        SAMPLE_AGG_STAT.out.aggregated_csv,
        SAMPLE_AGG_V.out.aggregated_csv
        )

    TCRDIST3_MATRIX(
        sample_map,
        params.matrix_sparsity,
        params.distance_metric,
        file(params.db_path)
    )

    TCRDIST3_MATRIX.out.max_matrix_value
        .map { tcrdist_xmax -> tcrdist_xmax.text.trim().toDouble() }
        .collect()
        .map { values -> values.max() }
        .set { global_x_max_value }
    global_x_max_value.view { global_xmax -> "Global x max matrix value: $global_xmax" }

    TCRDIST3_HISTOGRAM_CALC( 
        TCRDIST3_MATRIX.out.tcrdist_output,
        params.matrix_sparsity,
        params.distance_metric,
        global_x_max_value
    )

    TCRDIST3_HISTOGRAM_CALC.out.max_histogram_count
        .map { tcrdist_ymax -> tcrdist_ymax.text.trim().toDouble() }
        .collect()
        .map { values -> values.max() }
        .set { global_y_max_value }
    global_y_max_value.view { global_ymax -> "Global y max matrix value: $global_ymax" }

    TCRDIST3_HISTOGRAM_PLOT( 
        TCRDIST3_HISTOGRAM_CALC.out.histogram_data,
        global_y_max_value
    )

    cdr3_pgen_file = cdr3_pgen.first()
    OLGA_SAMPLE_MERGE ( sample_map,
        cdr3_pgen_file )

    OLGA_SAMPLE_MERGE.out.olga_xmin
        .map { xmin -> xmin.text.trim().toDouble() }
        .collect()
        .map { values -> values.min() }
        .set { olga_x_min_value }
    olga_x_min_value.view { olga_xmin -> "Olga x min matrix value: $olga_xmin" }

    OLGA_SAMPLE_MERGE.out.olga_xmax
        .map { xmax -> xmax.text.trim().toDouble() }
        .collect()
        .map { values -> values.max() }
        .set { olga_x_max_value }
    olga_x_max_value.view { olga_xmax -> "Olga x max matrix value: $olga_xmax" }

    OLGA_HISTOGRAM_CALC ( OLGA_SAMPLE_MERGE.out.olga_pgen, olga_x_min_value, olga_x_max_value )

    OLGA_HISTOGRAM_CALC.out.olga_ymax
        .map { ymax -> ymax.text.trim().toDouble() }
        .collect()
        .map { values -> values.max() }
        .set { olga_y_max_value }
    olga_y_max_value.view { olga_ymax -> "Olga y max matrix value: $olga_ymax" }

    OLGA_HISTOGRAM_PLOT( OLGA_HISTOGRAM_CALC.out.olga_histogram, olga_y_max_value )

    OLGA_WRITE_MAX(
        olga_x_min_value,
        olga_x_max_value,
        olga_y_max_value
    )

    CONVERGENCE ( sample_map )

    TCRPHENO ( sample_map )

    VDJDB_GET ()

    VDJDB_VDJMATCH (sample_map, VDJDB_GET.out.ref_db)

}