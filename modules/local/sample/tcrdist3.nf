process TCRDIST3_MATRIX {
    tag "${sample_meta.sample}"

    cpus {
            task.memory > 256.GB ? 16 * task.attempt :
            task.memory >  64.GB ?  8 * task.attempt :
            task.memory >   4.GB ?  4 * task.attempt :
                                    2 * task.attempt
    }


    memory {
        count_table.size() > 26 * 1024**2 ? 512.GB * task.attempt :
        count_table.size() > 20 * 1024**2 ? 256.GB * task.attempt :
        count_table.size() > 10 * 1024**2 ? 128.GB * task.attempt :
        count_table.size() >  4 * 1024**2 ?  64.GB * task.attempt :
        count_table.size() >  2 * 1024**2 ?  16.GB * task.attempt :
                                              4.GB * task.attempt
    }

    input:
    tuple val(sample_meta), path(count_table)
    val matrix_sparsity
    val distance_metric
    path ref_db

    output:
    tuple val(sample_meta), path("${sample_meta.sample}_distance_matrix.*"), emit: 'tcrdist_output'
    path "${sample_meta.sample}_clone_df.csv", emit: 'clone_df'
    path "matrix_maximum_value.txt", emit: 'max_matrix_value'

    script:
    """
    tcrdist3_matrix.py ${count_table} ${sample_meta.sample} ${matrix_sparsity} ${distance_metric} ${ref_db} ${task.cpus}
    """
}
