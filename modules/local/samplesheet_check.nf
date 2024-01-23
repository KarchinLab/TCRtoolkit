process SAMPLESHEET_CHECK {
    tag "${sample_table}"

    container "domebraccia/bulktcr:1.0-beta"

    publishDir "${params.output_dir}/input_check", mode: 'copy'

    input:
    path sample_table

    output:
    path 'sample_table_utf8.csv'    , emit: sample_utf8
    path 'sample_check.txt'

    script: 
    """
    #!/bin/bash
    
    iconv -t utf-8 $sample_table > sample_table_utf8.csv

    csvstat sample_table_utf8.csv > sample_check.txt
    """

    stub:
    """
    #!/bin/bash

    touch sample_check.txt
    """
}
