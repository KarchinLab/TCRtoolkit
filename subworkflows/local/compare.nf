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

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// include { CHECK_INPUT   } from '../../modules/local/check_input.nf'
include { CALC_COMPARE  } from '../../modules/local/calc_compare.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow COMPARE {

    println("Welcome to the BULK TCRSEQ pipeline! -- COMPARE ")

    /////// =================== CHECK INPUT ===================  ///////
    // CHECK_INPUT(
    //     file(params.sample_table), 
    //     file(params.patient_table) 
    //     )

    /////// =================== CALC COMPARE ==================  ///////
    // CALC_COMPARE( CHECK_INPUT.out.sample_utf8,
    //               CHECK_INPUT.out.patient_utf8 )


    /////// =================== PLOT COMPARE ===================  ///////
    
}