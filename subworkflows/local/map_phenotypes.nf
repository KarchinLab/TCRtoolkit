//
// Generate phenotype-specific samplesheet and transform phenotype files
//

include { GENERATE_PHENO_SAMPLESHEET } from '../../modules/local/samplesheet/generate_pheno_samplesheet.nf'
import groovy.json.JsonOutput

workflow MAP_PHENOTYPES {

    take:
    ch_phenotype_files_raw 
    original_samplesheet   

    main:

    // --- 1. Define transformation logic ---
    ch_phenotype_files_transformed = ch_phenotype_files_raw
    .transpose()
    .map { meta, file ->
        def original_sample_id = meta.sample
        def phenotype = file.name
                            .replaceFirst("^${original_sample_id}_", "")
                            .replaceFirst("_pseudobulk_phenotype.tsv\$", "")

        def new_sample_id = original_sample_id + '_' + phenotype

        def new_meta = meta.clone()
        new_meta.sample = new_sample_id
        new_meta.phenotype = phenotype
        new_meta.original_sample = original_sample_id

        new_meta.file = file.toString() 
        return [ new_meta, file ]               

    }

    // --- 2. Prepare inputs for samplesheet process ---
    ch_meta_json = ch_phenotype_files_transformed
        .map { meta, file -> meta }
        .collect() // Collect channel into a list
        .map { collected_list ->
            def json_string = JsonOutput.toJson(collected_list)
            def json_file = file("phenotype_meta.json") // Write to local temp file
            json_file.write(json_string)
            return json_file
        }

    // original_samplesheet_ = file(original_samplesheet)
    ch_script = file("$baseDir/bin/create_pheno_samplesheet.py")
    
    // --- 3. Call GENERATE_PHENO_SAMPLESHEET ---
    GENERATE_PHENO_SAMPLESHEET(
        original_samplesheet, // original_samplesheet_,
        ch_meta_json,
        ch_script
    )

    // View the output (for debugging)
    GENERATE_PHENO_SAMPLESHEET.out.samplesheet.view { file ->
        "--- START: Output of phenotype_samplesheet.csv ---\n" +
        file.text +
        "--- END: Output of phenotype_samplesheet.csv ---"
    }

    emit:
    files_transformed = ch_phenotype_files_transformed
    samplesheet_pheno = GENERATE_PHENO_SAMPLESHEET.out.samplesheet
}