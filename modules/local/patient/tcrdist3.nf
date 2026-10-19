/*
 * TCRDIST3_MATRIX
 *
 * One TCRdist3 distance matrix per patient, over the clones of all that patient's samples,
 * so similar TCRs are linked across timepoints and tissues (the scope GIANA and GLIPH2
 * cluster at). Each sample's block of the patient matrix is also written as
 * <sample>_distance_matrix.* / <sample>_clone_df.csv for the per-sample report sections.
 */
process TCRDIST3_MATRIX {
    tag "${patient}"

    cpus {
            task.memory > 256.GB ? 16 * task.attempt :
            task.memory >  64.GB ?  8 * task.attempt :
            task.memory >   4.GB ?  4 * task.attempt :
                                    2 * task.attempt
    }


    memory {
        def input_size = count_tables.sum { f -> f.size() }
        input_size > 26 * 1024**2 ? 512.GB * task.attempt :
        input_size > 20 * 1024**2 ? 256.GB * task.attempt :
        input_size > 10 * 1024**2 ? 128.GB * task.attempt :
        input_size >  4 * 1024**2 ?  64.GB * task.attempt :
        input_size >  2 * 1024**2 ?  16.GB * task.attempt :
                                       4.GB * task.attempt
    }

    input:
    tuple val(patient), val(samples), path(count_tables, arity: "1..*")
    val matrix_sparsity
    val distance_metric
    path ref_db

    output:
    tuple val(patient), path("${patient}_tcrdist_patient_matrix.*"), emit: 'patient_matrix'
    path "${patient}_tcrdist_patient_clones.csv", emit: 'patient_clones'
    path "*_distance_matrix.*", emit: 'sample_matrices'
    path "*_clone_df.csv", emit: 'sample_clone_dfs'
    path "matrix_maximum_value.txt", emit: 'max_matrix_value'

    script:
    def inputs = [samples, count_tables].transpose().collect { s, f -> "'${s}=${f}'" }.join(' ')
    """
    tcrdist3_matrix.py '${patient}' ${matrix_sparsity} ${distance_metric} ${ref_db} ${task.cpus} \\
        --input ${inputs}
    """
}
