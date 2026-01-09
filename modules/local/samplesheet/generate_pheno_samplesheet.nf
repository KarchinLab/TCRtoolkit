process GENERATE_PHENO_SAMPLESHEET {
    tag "Generating phenotype samplesheet"
    label 'process_single'
    // container "ghcr.io/karchinlab/tcrtoolkit:main"

    input:
    path original_samplesheet
    path meta_json 
    path py_script  // NOTE: Remove when in container

    output:
    path "phenotype_samplesheet.csv", emit: samplesheet

    script:
    def out_file = "phenotype_samplesheet.csv" // Define the output file variable

    """ 
    ${py_script} ${meta_json} ${out_file}
    """
}