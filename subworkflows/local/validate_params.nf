include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

workflow VALIDATE_PARAMS {   

    main:
    validateParameters()
    
    // Parameter summary logging is off: paramsSummaryLog fails under Nextflow 26.04.
    // log.info paramsSummaryLog(workflow)
}