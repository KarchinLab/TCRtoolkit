/*
 * MERGE_VDJ_OBJECT
 *
 * Builds a merged per-cell TCR object from the VDJ_QC contig tables, at both the pre- and
 * post-QC stage.
 *
 * Why this exists: on the VDJ-only route nothing emits a per-cell object. The pipeline
 * pools samples into pseudobulk, and BULK_TO_EXPORT then synthesizes a clonotype-level
 * export whose "cells" are reconstructed rather than real barcodes (4,372 synthetic rows
 * against 23,764 actual cells on the 8-sample test set). This process keeps real barcode
 * identity, so cross-sample clonal questions can be asked at the cell level.
 *
 * It runs on both routes. With a GEX object TCELL_INTEGRATION already produces the
 * authoritative Seurat, but the pre/post-QC pair here is still the cleanest way to see what
 * VDJ QC actually removed.
 *
 * Note the emitted Seurat carries a placeholder counts assay — cellranger vdj has no genes.
 * See bin/merge_vdj_object.R's header for what that does and does not permit.
 */
process MERGE_VDJ_OBJECT {
    tag "${project_name}"
    label 'process_medium'
    container "${params.sc_container}"

    publishDir "${params.outdir}/bridge/merged_vdj_object", mode: 'copy', overwrite: true

    input:
      path contigs_pre,   stageAs: 'contigs_before_qc.tsv'
      path contigs_post,  stageAs: 'contigs_after_qc.tsv'
      // Staged rather than called from ${projectDir}/bin so its contents join the task
      // hash; otherwise editing the script leaves -resume reusing the previous output.
      path merge_script,  stageAs: 'merge_vdj_object.R'
      // Post-merge GEX export when the full-SC route ran; NO_FILE otherwise. Used only to
      // flag which conserved receptors also reached the GEX analysis.
      path gex_export,    stageAs: 'in_gex_export.tsv'
      val  project_name

    output:
      path "pre_qc_*",  emit: pre_qc,  optional: true
      path "post_qc_*", emit: post_qc, optional: true
      path "*_cells.tsv", emit: cells, optional: true

    script:
    """
    # Each stage is independent: a failure on one must not cost the other, and a 0-byte
    # NO_FILE placeholder (either table absent) is skipped rather than treated as data.
    # stage:file pairs - the staged names are VDJ_QC's own (before/after), while the
    # output prefix is pre/post, so they must be mapped explicitly rather than derived.
    for pair in "pre:contigs_before_qc.tsv" "post:contigs_after_qc.tsv"; do
        stage="\${pair%%:*}"
        f="\${pair##*:}"
        if [ ! -s "\$f" ]; then
            echo "[merge_vdj] \$f missing or empty — skipping \${stage}-QC object"
            continue
        fi
        Rscript merge_vdj_object.R \\
            --contigs "\$f" \\
            --prefix "\${stage}_qc" \\
            --outdir . \\
            --sample-col "${params.vdj_meta_sample_col ?: 'sample'}" \\
            --gex-export in_gex_export.tsv \\
          || echo "[merge_vdj] WARN \${stage}-QC object failed; continuing"
    done
    """
}
