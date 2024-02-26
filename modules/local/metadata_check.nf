process METADATA_CHECK {
    tag "${patient_table}"

    container "domebraccia/bulktcr:1.0-beta"

    publishDir "${params.output_dir}/input_check", mode: 'copy'

    input:
    path patient_table

    output:
    path 'patient_table_utf8.csv'   , emit: patient_utf8
    path 'patient_check.txt'

    script: 
    """
    #!/bin/bash
    
    iconv -t utf-8 $patient_table > patient_table_utf8.csv

    csvstat patient_table_utf8.csv > patient_check.txt
    """

    stub:
    """
    #!/bin/bash

    touch patient_check.txt
    """
}
