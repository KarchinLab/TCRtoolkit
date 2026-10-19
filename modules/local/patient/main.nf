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
