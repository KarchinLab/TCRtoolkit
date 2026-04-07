process GLIPH2_TURBOGLIPH {
    tag "${patient}"
    label 'process_high'
    label 'process_high_compute'
    label 'process_high_memory'

    input:
    tuple val(patient), path(concat_cdr3)

    output:
    path "${patient}/all_motifs.txt", emit: 'all_motifs'
    path "${patient}/clone_network.txt", emit: 'clone_network'
    path "${patient}/cluster_member_details.txt", emit: 'cluster_member_details'
    path "${patient}/convergence_groups.txt", emit: 'convergence_groups'
    path "${patient}/global_similarities.txt", emit: 'global_similarities'
    path "${patient}/local_similarities.txt", emit: 'local_similarities'
    path "${patient}/parameter.txt", emit: 'gliph2_parameters'
    
    script:
    """
    Rscript - <<EOF
    #!/usr/bin/env Rscript

    library(turboGliph)

    # During testing, including TRBJ column was causing issues in clustering step. Removing and reinserting afterwards.
    df <- read.csv("$concat_cdr3", sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    df[,'patient'] <- df[,'sample']

    result <- turboGliph::gliph2(
        cdr3_sequences = df,
        result_folder = "./${patient}",
        lcminp = ${params.local_min_pvalue},
        sim_depth = ${params.simulation_depth},
        kmer_mindepth = ${params.kmer_min_depth},
        lcminove = ${params.local_min_OVE},
        all_aa_interchangeable = FALSE,
        n_cores = ${task.cpus}
    )

    df3 <- read.csv('${patient}/cluster_member_details.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
    df3[,'sample'] <- df3[,'patient']
    df3 <- merge(df3, df[, c("CDR3b", "TRBV", "sample", 'counts')], by = c("CDR3b", "TRBV", "sample", 'counts'), all.x = TRUE)
    df3 <- df3[, c('CDR3b', 'TRBV', 'TRBJ', 'counts', 'sample', 'tag', 'seq_ID', 'ultCDR3b')]
    write.table(df3, "${patient}/cluster_member_details.txt", sep = "\t", row.names = FALSE, quote = FALSE)
    EOF

    # Rename local_similarities file to standardize output name
    input_file="${patient}/local_similarities_*.txt"
    cat \$input_file > ${patient}/local_similarities.txt
    """
}
