process COMPARE_CALC {
    label 'process_single'
    container "ghcr.io/karchinlab/tcrtoolkit:main"
    
    input:
    path sample_utf8
    path all_sample_files

    output:
    path 'jaccard_mat.csv', emit: jaccard_mat
    path 'sorensen_mat.csv', emit: sorensen_mat
    path 'morisita_mat.csv', emit: morisita_mat

    script:
    """
    compare_calc.py \
        -s $sample_utf8 \
    """
}
