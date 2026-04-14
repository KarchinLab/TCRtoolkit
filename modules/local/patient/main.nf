process PATIENT_CONCATENATE {
    tag "${patient}"
    label "process_single"

    input:
    tuple val(patient), path(files)

    output:
    tuple val(patient), path("${patient}_cdr3.tsv"), emit: patient_cdr3

    script:
    """
    head -n 1 ${files[0]} > ${patient}_cdr3.tsv
    tail -n +2 -q ${files.join(' ')} >> ${patient}_cdr3.tsv
    """
}

process PATIENT_CALC {
    tag "${patient}"
    label 'process_single'

    input:
    tuple val(patient), path(concat_cdr3)

    output:
    path "${patient}/jaccard_mat.csv", emit: jaccard_mat
    path "${patient}/sorensen_mat.csv", emit: sorensen_mat
    path "${patient}/morisita_mat.csv", emit: morisita_mat

    script:
    """
    mkdir -p "${patient}"
    
    patient_calc.py \
        -i ${concat_cdr3} -p "${patient}" \
    """
}
