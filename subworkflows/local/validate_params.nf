include { validateParameters; paramsSummaryLog } from 'plugin/nf-schema'

workflow VALIDATE_PARAMS {   

    main:
    validateParameters()
    log.info paramsSummaryLog(workflow)
}