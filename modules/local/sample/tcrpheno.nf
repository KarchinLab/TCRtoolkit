process TCRPHENO {
    tag "${sample_meta.sample}"
    label 'process_low'
    container "ghcr.io/karchinlab/tcrtoolkit:main"

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

    df <- utils::read.csv("$count_table", sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
    df <- df %>%
        dplyr::filter(productive == 'true') %>% 
        dplyr::select(sequence_id, junction, junction_aa, v_call, j_call)

    df2 <- df %>%
        dplyr::rename(
            cell = sequence_id,
            TCRB_cdr3nt = junction,
            TCRB_cdr3aa = junction_aa,
            TCRB_vgene = v_call,
            TCRB_jgene = j_call
        )

    result <- tcrpheno::score_tcrs(df2,chain="b")
    result["sequence_id"] <- base::rownames(result)
    df3 <- dplyr::left_join(df, result, by = 'sequence_id')

    write.table(df3, "${sample_meta.sample}_tcrpheno.tsv", sep = "\t", row.names = FALSE, quote = FALSE, na = "")

    EOF
    """
}