process GIANA_CALC {
    label 'process_medium'

    input:
    path concat_cdr3
    val threshold
    val threshold_score
    val threshold_vgene

    output:
    path "VgeneScores.txt"
    path "giana_RotationEncodingBL62.txt"
    path "giana_EncodingMatrix.txt"
    path "giana.log"

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
        --Verbose \
        > giana.log 2>&1

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
    ' giana_RotationEncodingBL62.txt > tmp && mv tmp giana_RotationEncodingBL62.txt

    mv giana_RotationEncodingBL62.txt_EncodingMatrix.txt giana_EncodingMatrix.txt
    """
}
