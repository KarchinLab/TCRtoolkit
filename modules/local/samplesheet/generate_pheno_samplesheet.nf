process GENERATE_PHENO_SAMPLESHEET {
    label 'process_single'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    input:
    path(original_samplesheet)
    val(new_meta_list)

    output:
    path("phenotype_samplesheet.csv"), emit: samplesheet

    script:
    """
    // Read the original samplesheet and create a lookup map
    def lines = new File(original_samplesheet).readLines()
    def header_line = lines.first()
    def header_cols = header_line.split(',')

    // Find the 'sample' column index to be robust
    def sample_col_idx = header_cols.indexOf('sample')
    if (sample_col_idx == -1) {
        throw new Exception("Samplesheet header does not contain a 'sample' column.")
    }

    // Create a lookup map: 'Patient01_Base' -> 'Patient01_Base,Patient01,Base,tumor,...'
    def original_rows_map = lines.drop(1).collectEntries { line ->
        def cols = line.split(',')
        [ cols[sample_col_idx], line ]
    }

    // Create the new header (original + "phenotype")
    def new_header = header_line + ',phenotype'

    // Build the new rows from the input metadata list
    def new_rows = new_meta_list.collect { meta ->
        // Look up the original row string using the 'original_sample' key
        // e.g., 'Patient01_Base,Patient01,Base,tumor,...'
        def original_row_string = original_rows_map[meta.original_sample]

        // Split it: ['Patient01_Base', 'Patient01', 'Base', 'tumor', ...]
        def new_cols = original_row_string.split(',')

        // Update the sample ID column with the new one: 'Patient01_Base_CD4'
        new_cols[sample_col_idx] = meta.sample

        // Join back together and add the new phenotype value
        // 'Patient01_Base_CD4,Patient01,Base,tumor,...,CD4'
        return new_cols.join(',') + ',' + meta.phenotype
    }

    // Write the new samplesheet file
    def new_content = ([new_header] + new_rows).join('\\n')
    new File("phenotype_samplesheet.csv").write(new_content)
    """
}