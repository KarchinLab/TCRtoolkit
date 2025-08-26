process SAMPLESHEET_RESOLVE {
    label 'process_single'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    input:
    val(resolved_rows)     // List of tab-separated strings
    val(resolved_header)   // Comma-separated header line

    output:
    path "samplesheet_resolved.csv", emit: samplesheet_resolved

    script:
    """
    echo \"$resolved_header\" > samplesheet_resolved.csv

    for row in ${resolved_rows.collect{"\"${it}\""}.join(' ')}; do
        echo -e "\$row" >> samplesheet_resolved.csv
    done
    """
}