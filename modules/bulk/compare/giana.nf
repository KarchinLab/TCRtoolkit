process GIANA_CALC {
    tag "${patient}"
    label 'process_medium'
    publishDir "${params.outdir}/patient/giana", mode: 'copy'

    input:
    tuple val(patient), path(concat_cdr3)
    val threshold
    val threshold_score
    val threshold_vgene

    output:
    path "${patient}_VgeneScores.txt"
    path "${patient}_giana.txt"
    // path "giana_EncodingMatrix.txt"

    script:   
    """
    GIANA --file "${concat_cdr3}" \
        --output . \
        --outfile giana_RotationEncodingBL62.txt \
        --EncodingMatrix true \
        --threshold ${threshold} \
        --threshold_score ${threshold_score} \
        --threshold_vgene ${threshold_vgene} \
        --NumberOfThreads ${task.cpus} \
        --Verbose

    # Insert header after GIANA comments
    insert=\$(head -n 1 "${concat_cdr3}")
    insert=\$(echo "\$insert" | awk -F'\t' 'BEGIN{OFS="\t"} {
        out = \$1 OFS "cluster"
        for (i=2; i<=NF; i++) {
            out = out OFS \$i
        }
        print out
    }')

    awk -v insert="\$insert" '
    /^##/ { print; next }
    !inserted { print insert; inserted=1 }
    { print }
    ' giana_RotationEncodingBL62.txt > ${patient}_giana.txt

    # mv giana_RotationEncodingBL62.txt_EncodingMatrix.txt giana_EncodingMatrix.txt
    mv VgeneScores.txt ${patient}_VgeneScores.txt
    """
}
