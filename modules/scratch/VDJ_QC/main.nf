process VDJ_QC {

    tag "Running VDJ QC - ${project_name}"
    label 'process_medium'

    container "${params.container}"

    publishDir "${params.outdir}/VDJ_QC",mode: 'copy', overwrite: true

    input:
        path(notebook)
        path(sample_sheet)
        path (input_annotated_object)
        val(project_name)
        // Cirro/S3-backed datasets can't stage a whole outs/ directory - only
        // individual files (see .cirro_singlecell_*/preprocess.py). These are
        // empty lists for local/on-prem runs (params.contigs_files etc. unset),
        // in which case the stage_vdj_files.py step below is a no-op and
        // sample_sheet's own path column (a real local directory) is used as-is.
        path(contigs_files, stageAs: 'staged_vdj/contigs_input*/*')
        path(clonotypes_files, stageAs: 'staged_vdj/clonotypes_input*/*')
        path(metrics_files, stageAs: 'staged_vdj/metrics_input*/*')


    output:
        path("VDJ_QC_analysis.html"),               emit: report_html
        path("VDJ_QC/tables/contigs_after_qc.tsv"), emit: contigs_after_qc
        path("VDJ_QC/tables/*"),                    emit: tables
        path("VDJ_QC/figures/*"),                   emit: figures

    script:
        """
        mkdir -p VDJ_QC
        mkdir -p .cache/quarto
        export XDG_CACHE_HOME="\$PWD/.cache"
        export QUARTO_CACHE_DIR="\$PWD/.cache/quarto"
        export XDG_DATA_HOME="\$PWD/.cache"
        export QUARTO_PRINT_STACK=true
        export HOME="\$PWD"

        stage_vdj_files.py \
          --sample-sheet "${sample_sheet}" \
          --contigs-dir staged_vdj/contigs_input \
          --clonotypes-dir staged_vdj/clonotypes_input \
          --metrics-dir staged_vdj/metrics_input \
          --out-sample-sheet resolved_sample_sheet.csv

        quarto render ${notebook} \
          -P sample_sheet="resolved_sample_sheet.csv" \
          -P sample_sheet_sample_col="sample" \
          -P sample_sheet_path_col="path" \
          -P metadata_file="${params.metadata_file}" \
          -P input_annotated_object="${input_annotated_object}" \
          -P outdir="VDJ_QC" \
          -P chain_mode="${params.vdj_chain_mode}" \
          -P require_productive=${params.vdj_require_productive} \
          -P require_high_conf=${params.vdj_require_high_conf} \
          -P require_full_length=${params.vdj_require_full_length} \
          -P min_umis=${params.vdj_min_umis} \
          -P min_reads=${params.vdj_min_reads} \
          -P keep_one_per_chain=${params.vdj_keep_one_per_chain} \
          -P keep_paired_only=${params.vdj_keep_paired_only} \
          -P keep_dual_alpha=${params.vdj_keep_dual_alpha} \
          -P keep_dual_beta=${params.vdj_keep_dual_beta} \
          -P cdr3_aa_min=${params.vdj_cdr3_aa_min} \
          -P cdr3_aa_max=${params.vdj_cdr3_aa_max} \
          -P drop_stop_or_frames=${params.vdj_drop_stop_or_frames} \
          -P report_label="${project_name} VDJ QC"
        """
}