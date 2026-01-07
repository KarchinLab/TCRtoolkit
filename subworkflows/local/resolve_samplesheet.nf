//
// Check input samplesheet and get read channels
//

include { SAMPLESHEET_RESOLVE } from '../../modules/local/samplesheet/samplesheet_resolve'

workflow RESOLVE_SAMPLESHEET {
    take:
    samplesheet_utf8
    sample_map_final

    main:
    // Write resolved samplesheet with absolute paths
    sample_map_final
        .map { _meta, f -> f }
        .collect()
        .set { all_sample_files }
    
    sample_map_final
        .map { meta, f ->
            def row = (meta.values() + [f.getName()]).join(',')
            return row
        }
        .collect()
        .set { resolved_rows }

    samplesheet_utf8
        .splitCsv(header: true, sep: ',')
        .first()
        .map { row -> 
            def header = row.keySet().findAll { header_col -> header_col != 'file' } + ['file']
            return header.join(',')  // <-- convert to string
        }
        .set { resolved_header }

    SAMPLESHEET_RESOLVE(
            samplesheet_utf8,
            resolved_rows,
            resolved_header
        )
        .samplesheet_resolved
        .set { samplesheet_resolved }

    emit:
    samplesheet_resolved    // input to comparison analysis
    all_sample_files        // pass files to comparison tasks to be read by resolved samplesheet
    // versions = SAMPLESHEET_CHECK.out.versions // channel: [ versions.yml ]
}
