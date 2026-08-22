process SAMPLE_CALC {
    tag "${sample_meta.sample}"
    label 'process_single'
    publishDir enabled: false

    input:
    tuple val(sample_meta), path(count_table), path(pre_filter_stats)

    output:
    path "sample_stats_${sample_meta.sample}.csv"  , emit: sample_csv
    path "v_family_${sample_meta.sample}.csv"      , emit: v_family_csv
    path "d_family_${sample_meta.sample}.csv"      , emit: d_family_csv
    path "j_family_${sample_meta.sample}.csv"      , emit: j_family_csv

    script:
    """
    sample_calc.py -s '${sample_meta.sample}' -c ${count_table} -p ${pre_filter_stats}
    """
}

process SAMPLE_CALC_PIVOT {
    label 'process_single'
    publishDir "${params.outdir}/sample", mode: 'copy'

    input:
    path v_family_long
    path d_family_long
    path j_family_long

    output:
    path "v_family.csv", emit: v_family_wide
    path "d_family.csv", emit: d_family_wide
    path "j_family.csv", emit: j_family_wide

    script:
    """
    python - <<EOF
import pandas as pd
import re

def gene_sort_key(name):
    if m := re.search(r'(\\d+)\$', name):
        return (0, int(m.group(1)), '')
    if m := re.search(r'([A-Z]+)\$', name):
        return (1, 0, m.group(1))
    return (2, 0, name)

def pivot_family(input_file, output_file):
    df = pd.read_csv(input_file)
    wide = df.pivot_table(index='sample', columns='gene', values='count', fill_value=0).reset_index()
    gene_cols = [c for c in wide.columns if c != 'sample']
    wide = wide[['sample'] + sorted(gene_cols, key=gene_sort_key)]
    wide[gene_cols] = wide[gene_cols].astype(int)
    wide.to_csv(output_file, index=False)

pivot_family("${v_family_long}", "v_family.csv")
pivot_family("${d_family_long}", "d_family.csv")
pivot_family("${j_family_long}", "j_family.csv")
EOF
    """
}