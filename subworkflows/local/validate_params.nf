include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

workflow VALIDATE_PARAMS {   

    main:
    validateParameters()
    
    // TODO: Disabled due to error after updating to Nextflow 26.04.3. Re-enable when nf-core/schema catches up.
    // log.info paramsSummaryLog(workflow)
}