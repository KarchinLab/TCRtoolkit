process GLIPH2_TURBOGLIPH {
    label 'process_high'
    label 'process_high_compute'
    label 'process_high_memory'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

    input:
    path concat_cdr3

    output:
    path "all_motifs.txt", emit: 'all_motifs'
    path "clone_network.txt", emit: 'clone_network'
    path "cluster_member_details.txt", emit: 'cluster_member_details'
    path "convergence_groups.txt", emit: 'convergence_groups'
    path "global_similarities.txt", emit: 'global_similarities'
    path "local_similarities.txt", emit: 'local_similarities'
    path "parameter.txt", emit: 'gliph2_parameters'
    
    script:
    """
    # R script starts here
    cat > run_gliph2.R <<EOF
    #!/usr/bin/env Rscript

    library(turboGliph)

    # During testing, including TRBJ column was causing issues in clustering step. Removing and reinserting afterwards.
    df <- read.csv("$concat_cdr3", sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    # df2 <- subset(df, select = c('CDR3b', 'TRBV', 'patient', 'counts'))

    result <- turboGliph::gliph2(
        cdr3_sequences = df,
        result_folder = "./",
        lcminp = ${params.local_min_pvalue},
        sim_depth = ${params.simulation_depth},
        kmer_mindepth = ${params.kmer_min_depth},
        lcminove = ${params.local_min_OVE},
        all_aa_interchangeable = FALSE,
        n_cores = ${task.cpus}
    )

    df3 <- read.csv('cluster_member_details.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
    df3 <- merge(df3, df[, c("CDR3b", "TRBV", "patient", "TRBJ", 'counts')], by = c("CDR3b", "TRBV", "patient", 'counts'), all.x = TRUE)
    write.table(df3, "cluster_member_details.txt", sep = "\t", row.names = FALSE, quote = FALSE)

    EOF

    # Run the R script
    Rscript run_gliph2.R

    # Rename local_similarities file to standardize output name
    input_file="local_similarities_*.txt"
    cat \$input_file > local_similarities.txt
    """
}

process GLIPH2_PLOT {
    label 'process_low'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

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
        -P project_dir:$projectDir \
        -P results_dir:'./' \

        # -P clusters:$cluster_member_details \
        # -P cluster_stats:$convergence_groups \
        --to html
    """
}
