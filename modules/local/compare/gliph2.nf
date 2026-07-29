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
    // Patient-prefixed copy so collected files stay unique across patients (mirrors GIANA's
    // ${patient}_giana.txt). Consumed by the single-cell CLUSTER_TO_SC bridge.
    path "${patient}_cluster_member_details.txt", emit: 'cluster_member_details_named'

    script:
    """
    Rscript - <<EOF
    #!/usr/bin/env Rscript

    library(turboGliph)

    df <- read.csv("$concat_cdr3", sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    colnames(df)[colnames(df) == "junction_aa"]     <- "CDR3b"
    colnames(df)[colnames(df) == "v_call"]          <- "TRBV"
    colnames(df)[colnames(df) == "j_call"]          <- "TRBJ"
    colnames(df)[colnames(df) == "duplicate_count"] <- "counts"
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

    # Patient-prefixed top-level copy for unique staging when collected across patients.
    cp ${patient}/cluster_member_details.txt ${patient}_cluster_member_details.txt
    """
}

process GLIPH2_PLOT {
    label 'process_low'

    input:
    path gliph2_report_template
    path(motifs)
    path(clone_network)
    path(cluster_member_details)
    path(convergence_groups)
    path(global_similarities)
    path(local_similarities)
    path(parameter)

    output:
    path 'gliph2_report.html'

    script:   
    """
    ## copy quarto notebook to output directory
    cp $gliph2_report_template gliph2_report.qmd

    ## render qmd report to html
    quarto render gliph2_report.qmd \
        -P project_name:$params.project_name \
        -P workflow_cmd:'$workflow.commandLine' \
        -P results_dir:'./' \

        # -P clusters:$cluster_member_details \
        # -P cluster_stats:$convergence_groups \
        --to html
    """
}
