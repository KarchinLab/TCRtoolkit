
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CALC_SAMPLE } from '../../modules/local/calc_sample.nf'
include { PLOT_SAMPLE } from '../../modules/local/plot_sample.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SAMPLE {

    take:
    sample_map
    meta_data

    main:
    CALC_SAMPLE( sample_map,
                 meta_data )

    CALC_SAMPLE.out.sample_csv
        .collectFile(name: 'sample_stats.csv', sort: true, 
                     storeDir: params.output_dir)
        .map { file -> 
            ["bash", "-c", "echo 'sample_id,patient_id,timepoint,origin,treatment,response,clinical_data,num_clones,num_TCRs,simpson_index,simpson_index_corrected,clonality,num_in,num_out,num_stop,pct_prod,pct_out,pct_stop,pct_nonprod,cdr3_avg_len,num_convergent,ratio_convergent' \
            | cat - $file > temp && mv temp $file"].execute()
            return file 
        }
        .set { sample_stats_csv }

    CALC_SAMPLE.out.v_family_csv
        .collectFile(name: 'v_family.csv', sort: true,
                     storeDir: params.output_dir)
        .set { v_family_csv }

    CALC_SAMPLE.out.sample_meta
        .collectFile(name: 'sample_meta.csv', sort: true)
        .set { sample_meta_csv }
    
    /////// =================== PLOT SAMPLE ===================  ///////
    PLOT_SAMPLE(
        file(params.sample_table),
        sample_stats_csv,
        v_family_csv
        )
    
    // emit:
    // sample_stats_csv
    // v_family_csv
    // sample_meta_csv
    // versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}