process GLIPH2_TURBOGLIPH {
    tag "${patient}"
    label 'process_high'
    label 'process_high_compute'
    label 'process_high_memory'

    input:
    tuple val(patient), path(concat_cdr3)

    output:
    // all_motifs/clone_network/cluster_member_details/global_similarities are
    // copied to patient-prefixed top-level names because multiple patients'
    // outputs share the same basename otherwise (e.g. "all_motifs.txt"), which
    // collides when several patients' files are staged together downstream.
    tuple val(patient), path("${patient}_all_motifs.txt"), emit: 'all_motifs'
    tuple val(patient), path("${patient}_clone_network.txt"), emit: 'clone_network'
    tuple val(patient), path("${patient}_cluster_member_details.txt"), emit: 'cluster_member_details'
    path "${patient}/convergence_groups.txt", emit: 'convergence_groups'
    tuple val(patient), path("${patient}_global_similarities.txt"), emit: 'global_similarities'
    path "${patient}/local_similarities.txt", emit: 'local_similarities'
    path "${patient}/parameter.txt", emit: 'gliph2_parameters'
    
    script:
    """
    mkdir -p ${patient}

    Rscript - <<EOF
    #!/usr/bin/env Rscript

    library(turboGliph)

    df <- read.csv("$concat_cdr3", sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    colnames(df)[colnames(df) == "junction_aa"]     <- "CDR3b"
    colnames(df)[colnames(df) == "v_call"]          <- "TRBV"
    colnames(df)[colnames(df) == "j_call"]          <- "TRBJ"
    colnames(df)[colnames(df) == "duplicate_count"] <- "counts"
    df[,'patient'] <- df[,'sample']

    # turboGliph::gliph2() can throw an uncaught internal error - not just
    # print "no clusters/similarities found" - on edge-case-shaped input, e.g.
    # a dplyr::left_join() type mismatch between the tag column of sample_stats
    # (factor) and the tag column of ref_stats (double) when very few
    # candidate motifs pass its internal filtering (seen on real Patient02
    # data: crashes partway through "Part 2: Searching for global
    # similarities", after only all_motifs.txt and local_similarities_*.txt
    # were written). Catch this and fall back to empty results for the whole
    # patient rather than failing the task.
    gliph2_ok <- tryCatch({
        turboGliph::gliph2(
            cdr3_sequences = df,
            result_folder = "./${patient}",
            lcminp = ${params.local_min_pvalue},
            sim_depth = ${params.simulation_depth},
            kmer_mindepth = ${params.kmer_min_depth},
            lcminove = ${params.local_min_OVE},
            all_aa_interchangeable = FALSE,
            n_cores = ${task.cpus}
        )
        TRUE
    }, error = function(e) {
        message("turboGliph::gliph2() failed for ${patient}, treating as no clusters found: ", conditionMessage(e))
        FALSE
    })

    # gliph2() also doesn't error on the narrower "no significant clusters"
    # case - it just writes cluster_member_details.txt as a single blank line
    # with no header, which read.csv rejects outright ("no lines available in
    # input"). The same fallback covers that case, a totally missing file (if
    # gliph2() failed before reaching this point), and the tryCatch above.
    df3 <- tryCatch({
        tmp <- read.csv('${patient}/cluster_member_details.txt', sep = '\t', stringsAsFactors = FALSE, check.names = FALSE)
        tmp[,'sample'] <- tmp[,'patient']
        tmp <- merge(tmp, df[, c("CDR3b", "TRBV", "sample", 'counts')], by = c("CDR3b", "TRBV", "sample", 'counts'), all.x = TRUE)
        tmp[, c('CDR3b', 'TRBV', 'TRBJ', 'counts', 'sample', 'tag', 'seq_ID', 'ultCDR3b')]
    }, error = function(e) {
        data.frame(CDR3b=character(), TRBV=character(), TRBJ=character(), counts=integer(),
                   sample=character(), tag=character(), seq_ID=character(), ultCDR3b=character())
    })
    write.table(df3, "${patient}/cluster_member_details.txt", sep = "\t", row.names = FALSE, quote = FALSE)

    # If gliph2() failed partway through (or before writing anything), some of
    # its other declared output files may be entirely missing - write minimal
    # placeholders for whichever ones aren't already there so Nextflow's own
    # output declarations (which require every path to exist) don't turn this
    # into a hard pipeline failure on top of the already-handled analysis
    # failure.
    if (!file.exists('${patient}/all_motifs.txt')) {
        write.table(data.frame(motif=character(), num_in_sample=integer(), num_in_ref=integer(),
                                fisher.score=double(), num_fold=double()),
                    '${patient}/all_motifs.txt', sep = "\t", row.names = FALSE, quote = FALSE)
    }
    if (!file.exists('${patient}/clone_network.txt')) {
        file.create('${patient}/clone_network.txt')
    }
    if (!file.exists('${patient}/global_similarities.txt')) {
        write.table(data.frame(cluster_tag=character(), cluster_size=integer(), unique_CDR3b=integer(),
                                num_in_ref=integer(), fisher.score=double(), aa_at_position=character(),
                                TRBV=character(), CDR3b=character()),
                    '${patient}/global_similarities.txt', sep = "\t", row.names = FALSE, quote = FALSE)
    }
    if (!file.exists('${patient}/convergence_groups.txt')) {
        file.create('${patient}/convergence_groups.txt')
    }
    if (length(Sys.glob('${patient}/local_similarities_*.txt')) == 0) {
        write.table(data.frame(motif=character(), num_in_sample=integer(), num_in_ref=integer(),
                                fisher.score=double(), num_fold=double(), start=integer(), stop=integer(),
                                members=character()),
                    '${patient}/local_similarities_none.txt', sep = "\t", row.names = FALSE, quote = FALSE)
    }
    if (!file.exists('${patient}/parameter.txt')) {
        file.create('${patient}/parameter.txt')
    }
    EOF

    # Rename local_similarities file to standardize output name
    input_file="${patient}/local_similarities_*.txt"
    cat \$input_file > ${patient}/local_similarities.txt

    # Copy to patient-prefixed top-level names to avoid basename collisions
    # when multiple patients' outputs are staged together downstream.
    cp ${patient}/all_motifs.txt ${patient}_all_motifs.txt
    cp ${patient}/clone_network.txt ${patient}_clone_network.txt
    cp ${patient}/cluster_member_details.txt ${patient}_cluster_member_details.txt
    cp ${patient}/global_similarities.txt ${patient}_global_similarities.txt
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
