/*
 * SC_SAMPLE_STATS  (bridge)
 *
 * Synthesizes the pre-filter-stats sidecar that main's SAMPLE_CALC expects, for
 * single-cell-derived pseudobulk data. SC pseudobulk is productive-only by
 * construction (junction rebuilt from paired VDJ / GEX-annotated cells), so:
 *     total_clones = productive_clones = n rows ; nonproductive_clones = 0
 *
 * This lets the single-cell modality reuse main's unmodified 4-arg SAMPLE
 * subworkflow (Option A). Format matches ANNOTATE_PROCESS's pre_filter_stats:
 *     columns: sample,total_clones,productive_clones,nonproductive_clones (one row).
 */
process SC_SAMPLE_STATS {
    tag "${sample_meta.sample}"
    label 'process_low'
    publishDir enabled: false

    input:
    tuple val(sample_meta), path(count_table)

    output:
    tuple val(sample_meta), path("${sample_meta.sample}_pre_filter_stats.csv"), emit: "pre_filter_stats"

    script:
    """
    python - <<EOF
import pandas as pd

df = pd.read_csv("${count_table}", sep='\\t')
total = len(df)

pd.DataFrame({
    "sample": ["${sample_meta.sample}"],
    "total_clones": [total],
    "productive_clones": [total],
    "nonproductive_clones": [0],
}).to_csv("${sample_meta.sample}_pre_filter_stats.csv", index=False)
EOF
    """
}
