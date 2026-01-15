//
// Check input samplesheet and get read channels
//

workflow RESOLVE_SAMPLESHEET_PHENO {
    take:
    sample_map_final

    main:
    // Write resolved samplesheet with absolute paths
    sample_map_final
        .map { _meta, f -> f }
        .collect()
        .set { all_sample_files }

    emit:
    all_sample_files        // pass files to comparison tasks to be read by resolved samplesheet
}
