nextflow.enable.dsl = 2

process GLIPH2 {
    tag "${project_name}"
    label 'process_medium'
    container "syedsazaidi/scratch-gliph2:latest"

    publishDir "${params.outdir}/GLIPH2", mode: 'copy', overwrite: true

    input:
    path seurat_rds
    path export_cells
    path qmd
    path gliph_reference_bundle
    val  project_name
    

    output:
    path "GLIPH2_Report.html", emit: report_html
    path "GLIPH2_Report/data/seurat_with_GLIPH2.rds", emit: seurat_with_gliph2
    path "GLIPH2_Report/tables/gliph2_export_cells.tsv", emit: export_cells
    path "GLIPH2_Report/tables/*", emit: tables
    path "GLIPH2_Report/figures/*", emit: figures

    script:
    def gliph_clusters = params.gliph_clusters_file ?: ""
    def gliph_motifs   = params.gliph_motifs_file ?: ""
    def gliph_summary  = params.gliph_summary_file ?: ""

    """
    set -euo pipefail

    mkdir -p GLIPH2_Report

    test -f ${gliph_reference_bundle} || { echo "Missing GLIPH2 reference bundle: ${gliph_reference_bundle}" ; exit 1; }

    quarto render ${qmd} \
      -P seurat_rds="${seurat_rds}" \
      -P tcr_export_cells_file="${export_cells}" \
      -P gliph_clusters_file="${gliph_clusters}" \
      -P gliph_motifs_file="${gliph_motifs}" \
      -P gliph_summary_file="${gliph_summary}" \
      -P outdir="GLIPH2_Report" \
      -P gliph_ref_mode="${params.gliph_ref_mode}" \
      -P refdb_beta="${params.refdb_beta}" \
      -P gliph_reference_bundle="${gliph_reference_bundle}" \
      -P label_col="${params.label_col}" \
      -P sample_col="${params.sample_col}" \
      -P patient_col="${params.patient_col}" \
      -P condition_col="${params.condition_col}" \
      -P timepoint_col="${params.timepoint_col}" \
      -P batch_col="${params.batch_col}" \
      -P reduction_use="${params.reduction_use}" \
      -P make_umap_if_missing=${params.make_umap_if_missing} \
      -P umap_dims_max=${params.umap_dims_max} \
      -P umap_nfeatures=${params.umap_nfeatures} \
      -P raster_large_umap=${params.raster_large_umap} \
      -P clone_match_mode="${params.gliph_clone_match_mode}" \
      -P derive_clone_key_from_CTaa=${params.gliph_derive_clone_key_from_CTaa} \
      -P min_cells_per_group=${params.gliph_min_cells_per_group} \
      -P min_cluster_size_plot=${params.gliph_min_cluster_size_plot} \
      -P top_n_clusters=${params.gliph_top_n_clusters} \
      -P top_n_motifs=${params.gliph_top_n_motifs} \
      -P report_label="${project_name} GLIPH2"
    """
}


