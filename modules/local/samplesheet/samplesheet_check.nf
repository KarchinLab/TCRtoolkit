process SAMPLESHEET_CHECK {
    tag "${samplesheet}"
    label 'process_single'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    input:
    path samplesheet

    output:
    path 'samplesheet_utf8.csv'    , emit: samplesheet_utf8
    path 'samplesheet_stats.csv'

    script: 
    """
    samplesheet.py -s $samplesheet
    """

    stub:
    """
    #!/bin/bash

    touch samplesheet_utf8.csv
    touch samplesheet_stats.txt
    """
}
