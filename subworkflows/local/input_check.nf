//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_CHECK } from '../../modules/local/samplesheet_check'
include { METADATA_CHECK    } from '../../modules/local/metadata_check'

workflow INPUT_CHECK {
    take:
    samplesheet // file: /path/to/samplesheet.csv
    meta_data // file: /path/to/patientsheet.csv TODO: change to meta_data.csv

    main:

    SAMPLESHEET_CHECK( samplesheet )
        .sample_utf8
        .set { sample_utf8 }
    
    sample_utf8
        .splitCsv(header: true, sep: ',')
        .map { row -> 
            meta_map = [row.sample_id, row.patient_id, row.timepoint, row.origin] 
            [meta_map, file(row.file_path)]}
        .set { sample_map }
    
    view(sample_map)

    METADATA_CHECK( meta_data )
        .patient_utf8
        .set { meta_data }

    emit:
    sample_map
    meta_data
    sample_utf8
    // versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}
