include { SAMPLESHEET_CHECK } from '../../modules/bulk/samplesheet/samplesheet_check'

workflow INPUT_CHECK {
    take:
    samplesheet

    main:
    SAMPLESHEET_CHECK( samplesheet )
        .samplesheet_utf8
        .set { samplesheet_utf8 }

    samplesheet_utf8
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta     = row.findAll { k, _v -> k != 'file' }
            def file_obj = file(row.file)
            return [meta, file_obj]
        }
        .set { sample_map }

    emit:
    sample_map
    samplesheet_utf8
}
