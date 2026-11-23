#!/usr/bin/env Rscript
# SC_MARKER_EXPRESSION
#
# Mean expression of the canonical T-cell marker panel per cell state and clone group, from
# the annotated Seurat object. This is the only part of the single-cell reporting that needs
# the expression matrix: computing it once here lets the Master Summary (and anyone reading
# results/marker_expression.tsv) work from a few-KB table instead of loading the object.
#
#   sc_marker_expression.R --seurat_rds <rds> --export_cells <tsv> --out marker_expression.tsv
#
# Output, long format, one row per gene x state x group:
#   gene, annot, group (expanded: clone of >= 2 cells in its sample; singleton; no_tcr),
#   n_cells, mean_logexpr (mean log-normalised expression), frac_expressing (> 0)
# Means of several groups combine as sum(mean_logexpr * n_cells) / sum(n_cells).

suppressPackageStartupMessages({
    library(Seurat)
    library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
    i <- match(flag, args)
    if (is.na(i) || i == length(args)) default else args[[i + 1]]
}
rds_path    <- get_arg("--seurat_rds")
export_path <- get_arg("--export_cells")
out_path    <- get_arg("--out", "marker_expression.tsv")
genes_arg   <- get_arg("--genes", "")

genes <- if (nzchar(genes_arg)) strsplit(genes_arg, ",", fixed = TRUE)[[1]] else c(
    "CD8A", "CD8B", "CD4", "SELL", "CCR7", "TCF7", "IL7R", "GZMK", "GZMA", "GZMB",
    "PRF1", "GNLY", "NKG7", "IFNG", "CX3CR1", "KLRG1", "ITGAE", "CD69", "ZNF683",
    "PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "TOX", "CXCL13", "ENTPD1",
    "KLRB1", "SLC4A10", "ZBTB16", "MKI67", "TOP2A", "FOXP3", "TRDC")

empty_out <- function(why) {
    message("[sc_marker_expression] ", why, "; writing an empty table")
    write.table(data.frame(gene = character(), annot = character(), group = character(),
                           n_cells = integer(), mean_logexpr = numeric(), frac_expressing = numeric()),
                out_path, sep = "\t", quote = FALSE, row.names = FALSE)
    quit(save = "no", status = 0)
}

if (is.null(rds_path) || !file.exists(rds_path) || file.size(rds_path) == 0) empty_out("no Seurat object")
cells <- read.delim(export_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("cell_id", "sample") %in% names(cells))) empty_out("export lacks cell_id / sample")

seu <- readRDS(rds_path)
cat(sprintf("OBJECT_LOAD\tprocess=SC_MARKER_EXPRESSION\tfile=%s\tfile_mb=%.0f\tobject_mb=%.0f\tcells=%d\tfeatures=%d\tr_max_used_mb=%.0f",
            basename(rds_path), file.size(rds_path) / 2^20, as.numeric(utils::object.size(seu)) / 2^20,
            ncol(seu), nrow(seu), sum(gc()[, 6])), "\n", file = stderr())

genes <- intersect(genes, rownames(seu))
keep  <- intersect(as.character(cells$cell_id), colnames(seu))
if (length(genes) < 2 || length(keep) < 50) empty_out("too few marker genes or shared cells")

get_layer <- function(l) tryCatch(SeuratObject::LayerData(seu, layer = l),
                                  error = function(e) tryCatch(Seurat::GetAssayData(seu, slot = l),
                                                               error = function(e) NULL))
dat <- suppressWarnings(get_layer("data"))
if (!is.null(dat) && all(dim(dat) > 0)) {
    mat <- dat[genes, keep, drop = FALSE]
} else {
    # Counts only: LogNormalize the panel against each cell's full library size, as
    # NormalizeData would.
    cnt <- get_layer("counts")[, keep, drop = FALSE]
    mat <- log1p(Matrix::t(Matrix::t(cnt[genes, , drop = FALSE]) / Matrix::colSums(cnt)) * 1e4)
}
rm(seu, dat); invisible(gc())

# Clone group: expansion judged within the sample, as in the reports.
md <- cells[match(keep, cells$cell_id), , drop = FALSE]
has_tcr <- if ("has_tcr" %in% names(md)) toupper(as.character(md$has_tcr)) %in% c("TRUE", "T", "1") else !is.na(md$clone_id)
size <- ave(seq_len(nrow(md)), md$sample, md$clone_id, FUN = length)
grp  <- ifelse(!has_tcr | is.na(md$clone_id), "no_tcr", ifelse(size >= 2, "expanded", "singleton"))
annot <- if ("annot" %in% names(md)) as.character(md$annot) else "Unannotated"
annot[is.na(annot) | !nzchar(annot)] <- "Unannotated"

key <- paste(annot, grp, sep = "\r")
out <- do.call(rbind, lapply(split(seq_along(key), key), function(idx) {
    m <- mat[, idx, drop = FALSE]
    data.frame(gene = rownames(m), annot = annot[idx[1]], group = grp[idx[1]], n_cells = length(idx),
               mean_logexpr = round(Matrix::rowMeans(m), 5),
               frac_expressing = round(Matrix::rowMeans(m > 0), 5), stringsAsFactors = FALSE)
}))
out <- out[order(out$annot, out$group, match(out$gene, genes)), ]
write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("[sc_marker_expression] ", nrow(out), " rows: ", length(genes), " genes x ",
        length(unique(annot)), " states x clone groups")
