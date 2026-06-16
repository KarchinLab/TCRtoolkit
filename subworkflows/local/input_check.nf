//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet

    main:

    // 1. Run samplesheet_check
    SAMPLESHEET_CHECK( samplesheet )
        .samplesheet_utf8
        .set { samplesheet_utf8 }

    emit:
    sample_map = samplesheet_utf8
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = row.findAll { k, _v -> k != 'file' }  // everything except the file column
            def file_obj = file(row.file)
            return [meta, file_obj]
        }         //input to sample-level analysis
    samplesheet_utf8
    // versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}
