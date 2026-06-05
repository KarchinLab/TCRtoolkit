process TCRPHENO {
    tag "${sample_meta.sample}"
    label 'process_low'

    input:
    tuple val(sample_meta), path(count_table)

    output:
    tuple val(sample_meta), path("${sample_meta.sample}_tcrpheno.tsv"), emit: 'tcrpheno_output'

    script:
    """
    Rscript - <<EOF
#!/usr/bin/env Rscript

library(dplyr)
library(tcrpheno)

df <- utils::read.csv("${count_table}", sep = "\t", stringsAsFactors = FALSE, check.names = FALSE) %>%
    dplyr::select(sequence_id, junction, junction_aa, v_call, j_call)

df2 <- df %>%
    dplyr::rename(
        cell = sequence_id,
        TCRB_cdr3nt = junction,
        TCRB_cdr3aa = junction_aa,
        TCRB_vgene = v_call,
        TCRB_jgene = j_call
    ) %>%
    tcrpheno::score_tcrs(chain = "b")
df2["sequence_index"] <- base::as.integer(base::rownames(df2))

df3 <- df %>%
    dplyr::select(sequence_id) %>%
    dplyr::mutate(sequence_index = dplyr::row_number()) %>%
    dplyr::left_join(df2, by = 'sequence_index') %>%
    dplyr::select(!sequence_index)

write.table(df3, "${sample_meta.sample}_tcrpheno.tsv", sep = "\t", row.names = FALSE, quote = FALSE, na = "")
EOF
    """
}