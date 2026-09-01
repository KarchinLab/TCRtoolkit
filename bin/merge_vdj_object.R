#!/usr/bin/env Rscript
#
# merge_vdj_object.R
#
# Builds a merged per-cell TCR object from Cell Ranger VDJ contigs, for runs where no
# gene-expression matrix exists (the VDJ-only route, or any `cellranger vdj` output).
#
# IMPORTANT — what this can and cannot be:
#   `cellranger vdj` produces receptor sequences per cell barcode and NO genes. A Seurat
#   object is fundamentally a counts matrix (genes x cells) with metadata attached, so a
#   genuine expression-bearing Seurat cannot be built from VDJ data alone. This script
#   therefore emits two things:
#
#     1. <prefix>_combineTCR.rds   scRepertoire::combineTCR() output — the canonical
#                                  VDJ-only merged object. This is the honest artifact and
#                                  is exactly what you hand to combineExpression() later if
#                                  a GEX object turns up.
#
#     2. <prefix>_seurat.rds       A Seurat object carrying one cell per barcode with all
#                                  TCR fields in @meta.data, backed by a PLACEHOLDER assay
#                                  of zeros. Provided only for downstream tools that demand
#                                  the Seurat class. Do NOT normalize, run PCA/UMAP or
#                                  cluster on it — there is no expression to analyse.
#
#     3. <prefix>_cells.tsv        The same per-cell table as a flat TSV.
#
# Usage:
#   merge_vdj_object.R --contigs <contigs.tsv> --prefix <pre_qc|post_qc> [--outdir .]

suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
})

# Base-R argument parsing: optparse is not in the single-cell container.
.args <- commandArgs(trailingOnly = TRUE)
.get <- function(flag, default = NA_character_) {
    i <- match(flag, .args)
    if (!is.na(i) && length(.args) > i) .args[i + 1] else default
}
opt <- list(
    contigs    = .get("--contigs"),
    prefix     = .get("--prefix", "vdj"),
    outdir     = .get("--outdir", "."),
    sample_col = .get("--sample-col", "sample")
)
if (is.na(opt$contigs)) stop("--contigs is required")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) message("[merge_vdj] ", ...)

contigs <- as.data.frame(fread(opt$contigs, sep = "\t"))
msg(sprintf("read %s contigs, %d columns", format(nrow(contigs), big.mark = ","), ncol(contigs)))

scol <- if (opt$sample_col %in% names(contigs)) opt$sample_col else "sample"
stopifnot(scol %in% names(contigs), "barcode" %in% names(contigs))

# Cell Ranger barcodes repeat across samples, so the cell key must include the sample.
contigs$cell_id <- paste0(contigs[[scol]], "_", sub("-1$", "", contigs$barcode))

# ── 1. scRepertoire combineTCR ────────────────────────────────────────────────
# combineTCR() parses contig rows positionally and expects Cell Ranger's own column
# names, so feed it the raw frame split by sample rather than a renamed one.
tcr_obj <- NULL
if (requireNamespace("scRepertoire", quietly = TRUE)) {
    contig_list <- split(contigs, contigs[[scol]])
    samples     <- names(contig_list)
    tcr_obj <- tryCatch({
        scRepertoire::combineTCR(contig_list, samples = samples,
                                 removeNA = FALSE, removeMulti = FALSE, filterMulti = FALSE)
    }, error = function(e) {
        msg("WARN combineTCR() failed: ", conditionMessage(e)); NULL
    })
    if (!is.null(tcr_obj)) {
        saveRDS(tcr_obj, file.path(opt$outdir, paste0(opt$prefix, "_combineTCR.rds")))
        msg(sprintf("wrote %s_combineTCR.rds (%d samples)", opt$prefix, length(tcr_obj)))
    }
} else {
    msg("WARN scRepertoire unavailable — skipping combineTCR output")
}

# ── 2. per-cell table ─────────────────────────────────────────────────────────
# One row per cell. Alpha/beta chains are collapsed separately so a cell with several
# contigs of the same chain keeps the highest-UMI one, matching Cell Ranger convention.
pick <- function(df, want) {
    d <- df[df$chain %in% want, , drop = FALSE]
    if (!nrow(d)) return(NULL)
    if ("umis" %in% names(d)) d <- d[order(-as.numeric(d$umis)), , drop = FALSE]
    d[1, , drop = FALSE]
}

by_cell <- contigs %>% group_split(cell_id)
rows <- lapply(by_cell, function(d) {
    a <- pick(d, c("TRA", "TRG")); b <- pick(d, c("TRB", "TRD"))
    g <- function(x, col) if (!is.null(x) && col %in% names(x)) as.character(x[[col]]) else NA_character_
    data.frame(
        cell_id    = d$cell_id[1],
        sample     = d[[scol]][1],
        barcode    = sub("-1$", "", d$barcode[1]),
        n_contigs  = nrow(d),
        n_chains   = length(unique(d$chain)),
        cdr3a      = g(a, "cdr3"),    cdr3b      = g(b, "cdr3"),
        cdr3a_nt   = g(a, "cdr3_nt"), cdr3b_nt   = g(b, "cdr3_nt"),
        trav       = g(a, "v_gene"),  trbv       = g(b, "v_gene"),
        traj       = g(a, "j_gene"),  trbj       = g(b, "j_gene"),
        trac       = g(a, "c_gene"),  trbc       = g(b, "c_gene"),
        umis_a     = if (!is.null(a) && "umis" %in% names(a)) as.numeric(a$umis) else NA_real_,
        umis_b     = if (!is.null(b) && "umis" %in% names(b)) as.numeric(b$umis) else NA_real_,
        clonotype  = if ("raw_clonotype_id" %in% names(d)) as.character(d$raw_clonotype_id[1]) else NA_character_,
        multi_alpha = sum(d$chain %in% c("TRA", "TRG")) > 1,
        multi_beta  = sum(d$chain %in% c("TRB", "TRD")) > 1,
        stringsAsFactors = FALSE
    )
})
cells <- bind_rows(rows)

# scRepertoire's CTaa convention, so this object joins cleanly to the GEX route's output.
cells$CTaa <- ifelse(is.na(cells$cdr3a) & is.na(cells$cdr3b), NA_character_,
               paste0(ifelse(is.na(cells$cdr3a), "NA", paste0("A:", cells$cdr3a)), "|",
                      ifelse(is.na(cells$cdr3b), "NA", paste0("B:", cells$cdr3b))))
cells$has_tcr    <- !is.na(cells$cdr3a) | !is.na(cells$cdr3b)
cells$paired_tcr <- !is.na(cells$cdr3a) & !is.na(cells$cdr3b)

# Clone size = cells sharing a CTaa within a sample.
cells <- cells %>%
    group_by(sample, CTaa) %>%
    mutate(clone_size = ifelse(is.na(CTaa), NA_integer_, dplyr::n())) %>%
    ungroup() %>%
    as.data.frame()

fwrite(cells, file.path(opt$outdir, paste0(opt$prefix, "_cells.tsv")), sep = "\t")
msg(sprintf("wrote %s_cells.tsv (%s cells, %d samples, %s paired)",
            opt$prefix, format(nrow(cells), big.mark = ","),
            length(unique(cells$sample)), format(sum(cells$paired_tcr), big.mark = ",")))

# ── 3. Seurat object with a placeholder assay ─────────────────────────────────
if (requireNamespace("Seurat", quietly = TRUE)) {
    ok <- tryCatch({
        # A Seurat object requires a counts matrix. There is no expression data here, so
        # this is a single all-zero feature row purely to satisfy the class contract. The
        # information lives entirely in @meta.data.
        # Two rows, not one: Seurat 5's Assay5 rejects a single-row layer with
        # "Layers must be two-dimensional objects" because it drops to a vector.
        m <- Matrix::Matrix(0, nrow = 2, ncol = nrow(cells), sparse = TRUE,
                            dimnames = list(c("PLACEHOLDER-no-GEX-1", "PLACEHOLDER-no-GEX-2"),
                                            cells$cell_id))
        md <- cells; rownames(md) <- md$cell_id
        seu <- Seurat::CreateSeuratObject(counts = m, meta.data = md,
                                          project = opt$prefix, assay = "TCR",
                                          min.cells = 0, min.features = 0)
        seu@misc$provenance <- list(
            source        = normalizePath(opt$contigs),
            stage         = opt$prefix,
            note          = paste("VDJ-only: no gene expression. The counts assay is a",
                                  "zero placeholder. Do not normalize, run PCA/UMAP or",
                                  "cluster on this object."),
            n_cells       = nrow(cells),
            n_samples     = length(unique(cells$sample)),
            created       = as.character(Sys.time())
        )
        saveRDS(seu, file.path(opt$outdir, paste0(opt$prefix, "_seurat.rds")))
        msg(sprintf("wrote %s_seurat.rds (%s cells, placeholder assay)",
                    opt$prefix, format(ncol(seu), big.mark = ",")))
        TRUE
    }, error = function(e) { msg("WARN Seurat object failed: ", conditionMessage(e)); FALSE })
}

# ── 4. summary ────────────────────────────────────────────────────────────────
summ <- cells %>%
    group_by(sample) %>%
    summarise(cells = dplyr::n(),
              with_tcr = sum(has_tcr), paired = sum(paired_tcr),
              pct_paired = round(100 * sum(paired_tcr) / dplyr::n(), 2),
              unique_clonotypes = dplyr::n_distinct(CTaa[!is.na(CTaa)]),
              multi_chain = sum(multi_alpha | multi_beta),
              .groups = "drop")
fwrite(summ, file.path(opt$outdir, paste0(opt$prefix, "_summary.tsv")), sep = "\t")
print(as.data.frame(summ))
msg("done")
