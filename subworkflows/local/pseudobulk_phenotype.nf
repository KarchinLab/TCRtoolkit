//
// Generate phenotype-specific samplesheet and transform phenotype files
//

include { SAMPLESHEET_PHENO } from '../../modules/local/samplesheet/samplesheet_pheno'

include { SAMPLE_CALC as SAMPLE_CALC_PHENO } from '../../modules/local/sample/sample_calc'
include { SAMPLE_AGGREGATE as SAMPLE_AGG_STAT_PHENO } from '../../modules/local/sample/sample_aggregate' 
include { SAMPLE_AGGREGATE as SAMPLE_AGG_V_PHENO } from '../../modules/local/sample/sample_aggregate'
include { SAMPLE_AGGREGATE as SAMPLE_AGG_D_PHENO } from '../../modules/local/sample/sample_aggregate'
include { SAMPLE_AGGREGATE as SAMPLE_AGG_J_PHENO } from '../../modules/local/sample/sample_aggregate'
include { SAMPLE_PLOT as SAMPLE_PLOT_PHENO } from '../../modules/local/sample/sample_plot'

include { COMPARE_CALC as COMPARE_CALC_PHENO } from '../../modules/local/compare/compare_calc'
include { COMPARE_PLOT as COMPARE_PLOT_PHENO } from '../../modules/local/compare/compare_plot'
include { ANNOTATE_CONCATENATE as COMPARE_CONCATENATE_PHENO } from '../../modules/local/annotate'

workflow PSEUDOBULK_PHENOTYPE {
    take:
    ch_phenotype_files_raw 
    original_samplesheet  
    levels

    main:

    // --- 1. Define transformation logic ---
    ch_phenotype_files_transformed = ch_phenotype_files_raw
        .transpose()
        .map { meta, file ->
            def original_sample_id = meta.sample
            def phenotype = file.name
                                .replaceFirst("^${original_sample_id}_", "")
                                .replaceFirst("_pseudobulk_phenotype.tsv\$", "")

            def new_sample_id = original_sample_id + '_' + phenotype

            def new_meta = meta.clone()
            new_meta.sample = new_sample_id
            new_meta.phenotype = phenotype
            new_meta.original_sample = original_sample_id

            new_meta.file = file.toString() 
            return [ new_meta, file ]               
        }

    // --- 2. Prepare inputs for samplesheet process ---
    ch_meta_list = ch_phenotype_files_transformed
        .map { meta, _file -> meta }
        .collect() // Collect channel into a list

    // --- 3. Call SAMPLESHEET_PHENO ---
    SAMPLESHEET_PHENO(
        original_samplesheet,
        ch_meta_list
    )
    samplesheet_pheno = SAMPLESHEET_PHENO.out.samplesheet
    if (levels.contains('sample')) {
        SAMPLE_CALC_PHENO( ch_phenotype_files_transformed )
        SAMPLE_CALC_PHENO.out.sample_csv.collect().set { sample_csv_files }
        SAMPLE_CALC_PHENO.out.v_family_csv.collect().set { v_family_csv_files }
        SAMPLE_CALC_PHENO.out.d_family_csv.collect().set { d_family_csv_files }
        SAMPLE_CALC_PHENO.out.j_family_csv.collect().set { j_family_csv_files }

        SAMPLE_AGG_STAT_PHENO(sample_csv_files, "sample_stats.csv")
        SAMPLE_AGG_V_PHENO(v_family_csv_files, "v_family.csv")
        SAMPLE_AGG_D_PHENO(d_family_csv_files, "d_family.csv")
        SAMPLE_AGG_J_PHENO(j_family_csv_files, "j_family.csv")

        SAMPLE_PLOT_PHENO(
            samplesheet_pheno,
            file(params.sample_stats_template),
            SAMPLE_AGG_STAT_PHENO.out.aggregated_csv,
            SAMPLE_AGG_V_PHENO.out.aggregated_csv
            )
    }

    ch_phenotype_files_transformed
        .map { _meta, f -> f }
        .collect()
        .set { all_sample_files }
    if (levels.contains('compare')) {
        COMPARE_CALC_PHENO( samplesheet_pheno,
            all_sample_files )

        COMPARE_CONCATENATE_PHENO( samplesheet_pheno,
            all_sample_files )

        COMPARE_PLOT_PHENO( samplesheet_pheno,
            COMPARE_CALC_PHENO.out.jaccard_mat,
            COMPARE_CALC_PHENO.out.sorensen_mat,
            COMPARE_CALC_PHENO.out.morisita_mat,
            file(params.compare_stats_template),
            params.project_name,
            all_sample_files
        )
    }

    emit:
    files_transformed = ch_phenotype_files_transformed
    samplesheet_pheno
}