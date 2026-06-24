include { SAMPLE_CALC }                                      from '../../modules/bulk/sample/sample_calc'
include { SAMPLE_PLOT }                                      from '../../modules/bulk/sample/sample_plot'
include { SAMPLE_AGGREGATE as SAMPLE_AGG_STAT }              from '../../modules/bulk/sample/sample_aggregate'
include { SAMPLE_AGGREGATE as SAMPLE_AGG_V }                 from '../../modules/bulk/sample/sample_aggregate'
include { SAMPLE_AGGREGATE as SAMPLE_AGG_D }                 from '../../modules/bulk/sample/sample_aggregate'
include { SAMPLE_AGGREGATE as SAMPLE_AGG_J }                 from '../../modules/bulk/sample/sample_aggregate'
include { TCRDIST3_MATRIX; TCRDIST3_HISTOGRAM_CALC; TCRDIST3_HISTOGRAM_PLOT } from '../../modules/bulk/sample/tcrdist3'
include { OLGA_SAMPLE_MERGE; OLGA_HISTOGRAM_CALC; OLGA_HISTOGRAM_PLOT }       from '../../modules/bulk/olga/main'
include { CONVERGENCE }                                      from '../../modules/bulk/sample/convergence'
include { TCRPHENO }                                         from '../../modules/bulk/sample/tcrpheno'
include { VDJDB_GET; VDJDB_VDJMATCH }                        from '../../modules/bulk/sample/tcrspecificity'

workflow SAMPLE {
    take:
    sample_map
    cdr3_pgen
    olga_stats

    main:
    SAMPLE_CALC( sample_map )

    SAMPLE_CALC.out.sample_csv.collect().set   { sample_csv_files }
    SAMPLE_CALC.out.v_family_csv.collect().set { v_family_csv_files }
    SAMPLE_CALC.out.d_family_csv.collect().set { d_family_csv_files }
    SAMPLE_CALC.out.j_family_csv.collect().set { j_family_csv_files }

    SAMPLE_AGG_STAT( sample_csv_files,   "sample_stats.csv" )
    SAMPLE_AGG_V   ( v_family_csv_files, "v_family.csv"     )
    SAMPLE_AGG_D   ( d_family_csv_files, "d_family.csv"     )
    SAMPLE_AGG_J   ( j_family_csv_files, "j_family.csv"     )

    SAMPLE_PLOT(
        file(params.samplesheet ?: params.sample_sheet),
        file(params.sample_stats_template),
        SAMPLE_AGG_STAT.out.aggregated_csv,
        SAMPLE_AGG_V.out.aggregated_csv
    )

    TCRDIST3_MATRIX( sample_map, params.matrix_sparsity, params.distance_metric, file(params.db_path) )

    TCRDIST3_MATRIX.out.max_matrix_value
        .map { f -> f.text.trim().toDouble() }
        .collect()
        .map { vals -> vals.max() }
        .set { global_x_max_value }

    TCRDIST3_HISTOGRAM_CALC( TCRDIST3_MATRIX.out.tcrdist_output, params.matrix_sparsity, params.distance_metric, global_x_max_value )

    TCRDIST3_HISTOGRAM_CALC.out.max_histogram_count
        .map { f -> f.text.trim().toDouble() }
        .collect()
        .map { vals -> vals.max() }
        .set { global_y_max_value }

    TCRDIST3_HISTOGRAM_PLOT( TCRDIST3_HISTOGRAM_CALC.out.histogram_data, global_y_max_value )

    OLGA_SAMPLE_MERGE( sample_map, cdr3_pgen.first(), olga_stats )
    OLGA_HISTOGRAM_CALC( OLGA_SAMPLE_MERGE.out.olga_pgen, olga_stats )

    OLGA_HISTOGRAM_CALC.out.olga_ymax
        .map { f -> f.text.trim().toDouble() }
        .collect()
        .map { vals -> vals.max() }
        .set { olga_y_max_value }

    OLGA_HISTOGRAM_PLOT( OLGA_HISTOGRAM_CALC.out.olga_histogram, olga_y_max_value )

    CONVERGENCE( sample_map )
    TCRPHENO( sample_map )

    VDJDB_GET()
    VDJDB_VDJMATCH( sample_map, VDJDB_GET.out.ref_db )
}
