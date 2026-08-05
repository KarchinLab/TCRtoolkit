/*
 * PSEUDOBULK_QC
 *
 * TCRtoolkit QC gate for single-cell → bulk (pseudobulk) TCR data.
 *
 * PSEUDOBULK_QC_CALC      computes per-sample n_clones / n_cells and a
 *                         PASS/FAIL flag against the configured thresholds,
 *                         plus non-redundant V gene-family usage.
 *
 * PSEUDOBULK_QC_AGGREGATE concatenates per-sample CSVs into a single table
 *                         (QC summary or V-family usage) for reporting.
 *
 * Diversity metrics (Shannon/Simpson/etc.) are deliberately NOT computed here —
 * the REPERTOIRE module already reports them on the same cell counts.
 */

process PSEUDOBULK_QC_CALC {
    tag "${meta.sample}"
    label 'process_single'
    container "${params.container}"
    publishDir enabled: false

    input:
    tuple val(meta), path(sample_file)
    val   min_clones
    val   min_cells

    output:
    tuple val(meta), path(sample_file), env('QC_PASS'), env('N_CLONES'), env('N_CELLS'), emit: scored
    path "qc_${meta.sample}.csv",      emit: qc_csv
    path "vfamily_${meta.sample}.csv", emit: v_family_csv

    script:
    """
    pseudobulk_qc.py ${sample_file} '${meta.sample}' ${min_clones} ${min_cells}

    QC_PASS=\$(awk -F, 'NR==2{print \$6}' 'qc_${meta.sample}.csv')
    N_CLONES=\$(awk -F, 'NR==2{print \$2}' 'qc_${meta.sample}.csv')
    N_CELLS=\$(awk -F, 'NR==2{print \$3}' 'qc_${meta.sample}.csv')
    """
}

process PSEUDOBULK_QC_AGGREGATE {
    tag "${output_file}"
    label 'process_low'
    container "${params.container}"
    publishDir "${params.outdir}/pseudobulk_qc", mode: params.publish_dir_mode, overwrite: true

    input:
    path csv_files
    val  output_file

    output:
    path output_file, emit: aggregated_csv

    script:
    """
    cat > aggregate.py <<EOF
    import sys
    import pandas as pd

    input_files = sys.argv[1:]
    merged = pd.concat([pd.read_csv(f) for f in input_files], axis=0, ignore_index=True)
    if 'sample' in merged.columns:
        merged = merged.sort_values('sample')
    merged.to_csv("${output_file}", index=False)
    EOF

    python3 aggregate.py ${csv_files}
    """
}
