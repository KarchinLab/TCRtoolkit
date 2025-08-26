process GIANA_CALC {
    label 'process_medium'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

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
    python3 - <<EOF
    input_file = "giana_RotationEncodingBL62.txt"
    concat_header_file = "${concat_cdr3}"

    with open(concat_header_file, 'r', encoding='utf-8') as f:
        header = f.readline().strip().split('\\t')
    header.insert(1, "cluster")
    header_line = '\\t'.join(header)

    with open(input_file, 'r', encoding='utf-8') as infile:
        lines = infile.readlines()

    with open(input_file, 'w', encoding='utf-8') as outfile:
        inserted = False
        for line in lines:
            if line.startswith("##"):
                outfile.write(line)
            elif not inserted:
                outfile.write(header_line + '\\n')
                outfile.write(line)
                inserted = True
            else:
                outfile.write(line)
    EOF

    mv giana_RotationEncodingBL62.txt_EncodingMatrix.txt giana_EncodingMatrix.txt
    """
}
