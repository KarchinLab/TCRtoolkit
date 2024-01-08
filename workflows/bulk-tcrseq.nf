#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    btc/bulktcrseq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/KarchinLab/bulk-tcrseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Validate pipeline parameters
def checkPathParamList = [ params.sample_table, params.patient_table ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.project_name) { project_name = params.project_name } else { exit 1, 'Project name not specified. Please, provide a --project_name=project_name !' }
if (params.sample_table) { sample_table = file(params.sample_table) } else { exit 1, 'Sample table not specified. Please, provide a --sample_table=/path/to/sample_table.csv !' }
if (params.patient_table) { patient_table = file(params.patient_table) } else { exit 1, 'Patient table not specified. Please, provide a --patient_table=/path/to/patient_table.csv !' }
if (params.output_dir) { output_dir = params.output_dir } else { exit 1, 'Output directory not specified. Please, provide a --output_dir=/path/to/output_dir !' }

log.info """\
BTC - BULK TCRSEQ PIPELINE
===================================
Project parameters:
- Project Name              : ${params.project_name}
- Sample Table (CSV)        : ${params.sample_table}
- Patient Meta-data (CSV)   : ${params.patient_table}
- Output Directory          : ${params.output_dir}
"""

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CHECK_INPUT } from '../modules/local/check_input.nf'
include { CALC_SAMPLE } from '../modules/local/calc_sample.nf'
include { PLOT_SAMPLE } from '../modules/local/plot_sample.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BULKTCR {

    println("Welcome to the BULK TCRSEQ pipeline!")

    /////// =================== CHECK INPUT ===================  ///////
    CHECK_INPUT(
        file(params.sample_table), 
        file(params.patient_table) 
        )

    CHECK_INPUT.out.sample_utf8
        .splitCsv(header: true, sep: ',')
        .map { row -> 
            meta_map = [row.sample_id, row.patient_id, row.timepoint, row.origin] 
            [meta_map, file(row.file_path)]}
        .set { sample_map }

    /////// =================== CALC SAMPLE ===================  ///////
    CALC_SAMPLE( sample_map )

    CALC_SAMPLE.out.sample_csv
        .collectFile(name: 'sample_stats.csv', sort: true, 
                     storeDir: params.output_dir)
        .set { sample_stats_csv }

    CALC_SAMPLE.out.v_family_csv
        .collectFile(name: 'v_family.csv', sort: true)
        .set { v_family_csv }

    CALC_SAMPLE.out.sample_meta
        // .view()
        .collectFile(name: 'sample_meta.csv', sort: true)
        .set { sample_meta_csv }
    
    /////// =================== PLOT SAMPLE ===================  ///////
    PLOT_SAMPLE(
        file(params.sample_table),
        sample_stats_csv,
        v_family_csv
        )
    
}