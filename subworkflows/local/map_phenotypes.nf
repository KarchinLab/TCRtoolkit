/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENERATE_PHENO_SAMPLESHEET } from '../../modules/local/samplesheet/generate_pheno_samplesheet'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN SUBWORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
 * subworkflows/local/map_phenotypes.nf
 *
 * This subworkflow:
 * 1. Takes the raw phenotype files
 * 2. Transforms them to create new metadata (streaming)
 * 3. Collects the new metadata to generate a new samplesheet (collecting)
 * 4. Emits both the transformed file channel and the new samplesheet
 */

workflow MAP_PHENOTYPES {
    take:
    ch_phenotype_files_raw 
    original_samplesheet 

    main:
    transformed_files_ch = ch_phenotype_files_raw.map { meta, file ->
        def original_sample_id = meta.sample
        def phenotype = file.name.removePrefix(original_sample_id + '_').removeSuffix('_pseudobulk_phenotype.tsv')
        def new_sample_id = original_sample_id + '_' + phenotype

        def new_meta = meta.clone()
        new_meta.sample = new_sample_id
        new_meta.phenotype = phenotype
        new_meta.original_sample = original_sample_id

        return [ new_meta, file ]
    }
    .into { transformed_for_emit; transformed_for_collect }

    // Collect all new metadata into a single list
    def ch_new_meta_list = transformed_for_collect.map { it[0] }.collect()

    // Call your module to generate the samplesheet (Collecting)
    GENERATE_PHENO_SAMPLESHEET(
        original_samplesheet,
        ch_new_meta_list
    )
    def ch_phenotype_samplesheet = GENERATE_PHENO_SAMPLESHEET.out.samplesheet

    emit:
    files_transformed = transformed_for_emit // Emits the streaming channel
    samplesheet_pheno = ch_phenotype_samplesheet     // Emits the single file
}