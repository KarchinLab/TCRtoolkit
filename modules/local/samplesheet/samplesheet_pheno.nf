process SAMPLESHEET_PHENO {
    tag "Generating phenotype samplesheet"
    label 'process_single'

    input:
    path original_samplesheet
    val meta_list 

    output:
    path "samplesheet_phenotype.csv", emit: samplesheet

    script:
    def json_meta = groovy.json.JsonOutput.toJson(meta_list)
    """
    cat << 'EOF' > phenotype_meta.json
    ${json_meta}
    EOF
    create_pheno_samplesheet.py phenotype_meta.json samplesheet_phenotype.csv
    """
}