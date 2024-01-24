process CALC_COMPARE {
    
    // tag
    // label
    container "domebraccia/bulktcr:1.0-beta"

    input:
    path sample_utf8
    path patient_utf8

    output:
    path 'jaccard_amat.csv', emit: jaccard_amat
    path 'sorensen_amat.csv', emit: sorensen_amat

    script:
    """
    python $projectDir/bin/calc_compare.py \
        -s $sample_utf8 \
        -p $patient_utf8
    """

}