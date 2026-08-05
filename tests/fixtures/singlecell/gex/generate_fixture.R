#!/usr/bin/env Rscript
# Generates a small but real, paired GEX + VDJ single-cell fixture for testing
# TCELL_INTEGRATION (and, downstream, CONGA / CONSENSUS_CLUSTERING / CLUSTER_TO_SC):
#   - seurat_annotated.rds : a real Seurat object (RNA counts + PCA/UMAP + metadata)
#   - contigs_after_qc.tsv : matching VDJ contigs (same barcodes/samples), in the
#     schema TCELL_INTEGRATION reads directly (sample, barcode, chain, cdr3,
#     v_gene, j_gene, raw_clonotype_id).
#
# Barcodes are shared between the two files so Seurat/scRepertoire's barcode
# harmonization can match them. Not biologically meaningful data - just
# numerically valid input real Seurat/scRepertoire code can run against.
#
# IMPORTANT: barcodes must be real Cell-Ranger-style 16nt-ACGT + "-1" strings,
# not human-readable IDs - scRepertoire::combineTCR() silently drops any cell
# whose barcode isn't in that format (confirmed by direct testing: a 4-row
# input with "PatientA_Base_CELL001-1"-style barcodes produced 0 output rows
# per sample; the same input with "AAACCTGAGAAACCAT-1"-style barcodes worked).
# combineTCR() also prepends "{sample}_" to its own internal barcode column
# regardless of the `ID` argument - Seurat's raw colnames should NOT have that
# prefix baked in, since TCELL_INTEGRATION's harmonize_barcodes() step exists
# specifically to reconcile that "{sample}_{barcode}" vs "{barcode}" mismatch.

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

set.seed(42)

out_dir <- dirname(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
if (length(out_dir) == 0 || out_dir == "") out_dir <- "."

samples <- c("PatientA_Base", "PatientA_Post", "PatientB_Base")
patient_of <- c(PatientA_Base = "PatientA", PatientA_Post = "PatientA", PatientB_Base = "PatientB")
timepoint_of <- c(PatientA_Base = "Base", PatientA_Post = "Post", PatientB_Base = "Base")
cells_per_sample <- 60

# CoNGA refuses to run below a minimum clonotype count ("ERROR too few
# clonotypes", confirmed via its own real log output) - it needs enough
# distinct clones to find meaningful GEX-TCR statistical correlations, not
# just enough cells. 20 distinct clones comfortably clears that.
n_clones <- 20
trbv_pool <- c("TRBV7-2","TRBV5-1","TRBV6-5","TRBV4-1","TRBV27","TRBV19","TRBV20-1",
               "TRBV6-1","TRBV9","TRBV28","TRBV11-2","TRBV12-3","TRBV14","TRBV15",
               "TRBV18","TRBV24-1","TRBV25-1","TRBV29-1","TRBV30","TRBV10-3")
trbj_pool <- c("TRBJ2-7","TRBJ1-2","TRBJ2-3","TRBJ1-3","TRBJ1-4","TRBJ2-1","TRBJ2-5",
               "TRBJ1-1","TRBJ1-5","TRBJ2-6")
trav_pool <- c("TRAV3","TRAV13-1","TRAV21","TRAV12-3","TRAV27","TRAV1-2","TRAV17",
               "TRAV19","TRAV38-1","TRAV26-1","TRAV9-2","TRAV14/DV4","TRAV35",
               "TRAV8-1","TRAV24","TRAV5","TRAV20","TRAV22","TRAV30","TRAV41")
traj_pool <- c("TRAJ13","TRAJ6","TRAJ48","TRAJ33","TRAJ23","TRAJ40","TRAJ22",
               "TRAJ34","TRAJ42","TRAJ54")

aa <- c("A","S","G","P","T","N","D","E","Q","K","R","L","V","I","F","Y","W","H")
random_cdr3 <- function(prefix, min_len = 11, max_len = 16) {
  mid <- paste(sample(aa, sample(min_len:max_len, 1), replace = TRUE), collapse = "")
  paste0(prefix, mid, "F")
}

beta_clones <- lapply(seq_len(n_clones), function(i) {
  c(random_cdr3("CASS"), trbv_pool[[((i - 1) %% length(trbv_pool)) + 1]], trbj_pool[[((i - 1) %% length(trbj_pool)) + 1]])
})
alpha_clones <- lapply(seq_len(n_clones), function(i) {
  c(random_cdr3("CA"), trav_pool[[((i - 1) %% length(trav_pool)) + 1]], traj_pool[[((i - 1) %% length(traj_pool)) + 1]])
})
stopifnot(length(unique(sapply(beta_clones, `[[`, 1))) == n_clones)

random_nt_barcode <- function() paste0(paste(sample(c("A","C","G","T"), 16, replace = TRUE), collapse = ""), "-1")

all_barcodes <- c()
all_meta <- data.frame()
contig_rows <- list()

for (s in samples) {
  # Weight so a handful of clones are visibly expanded, not perfectly uniform.
  clone_prob <- rev(seq_len(n_clones))
  clone_idx <- sample(seq_len(n_clones), cells_per_sample, replace = TRUE,
                       prob = clone_prob / sum(clone_prob))
  barcodes <- character(0)
  while (length(barcodes) < cells_per_sample) {
    cand <- random_nt_barcode()
    if (!(cand %in% all_barcodes) && !(cand %in% barcodes)) barcodes <- c(barcodes, cand)
  }
  all_barcodes <- c(all_barcodes, barcodes)

  # CONGA (unlike TCELL_INTEGRATION) has no fallback when no cell-type
  # annotation column is present - it hard-requires one of
  # predicted_labels/celltype/Annotation and errors otherwise
  # ("Could not resolve annotation label column."), since correlating GEX
  # cell state with TCR clusters is the whole point of the tool. Real usage
  # is expected to supply an already-annotated GEX object; synthesize a
  # plausible-enough placeholder here so the module can actually run.
  celltype <- sample(c("CD8 T cell", "CD4 T cell"), cells_per_sample, replace = TRUE, prob = c(0.6, 0.4))

  all_meta <- rbind(all_meta, data.frame(
    barcode = barcodes,
    orig.ident = s,
    patient_id = patient_of[[s]],
    condition = "tumor",
    timepoint = timepoint_of[[s]],
    celltype = celltype,
    stringsAsFactors = FALSE
  ))

  for (i in seq_len(cells_per_sample)) {
    ci <- clone_idx[i]
    clonotype_id <- paste0("clonotype", ci)
    beta_nt <- paste0("TGT", paste(rep("GCC", nchar(beta_clones[[ci]][1]) - 2), collapse = ""), "TTT")
    alpha_nt <- paste0("TGT", paste(rep("GCC", nchar(alpha_clones[[ci]][1]) - 2), collapse = ""), "TTT")
    contig_rows[[length(contig_rows) + 1]] <- data.frame(
      sample = s, barcode = barcodes[i], is_cell = "True",
      contig_id = paste0(barcodes[i], "_contig_1"), high_confidence = "True",
      length = 550, chain = "TRB",
      cdr3 = beta_clones[[ci]][1], cdr3_nt = beta_nt,
      v_gene = beta_clones[[ci]][2], j_gene = beta_clones[[ci]][3],
      d_gene = "TRBD1", c_gene = "TRBC1", full_length = "True",
      productive = "True", reads = 1200, umis = 4,
      raw_clonotype_id = clonotype_id,
      raw_consensus_id = paste0(clonotype_id, "_consensus_1"), stringsAsFactors = FALSE
    )
    contig_rows[[length(contig_rows) + 1]] <- data.frame(
      sample = s, barcode = barcodes[i], is_cell = "True",
      contig_id = paste0(barcodes[i], "_contig_2"), high_confidence = "True",
      length = 500, chain = "TRA",
      cdr3 = alpha_clones[[ci]][1], cdr3_nt = alpha_nt,
      v_gene = alpha_clones[[ci]][2], j_gene = alpha_clones[[ci]][3],
      d_gene = "None", c_gene = "TRAC", full_length = "True",
      productive = "True", reads = 1100, umis = 3,
      raw_clonotype_id = clonotype_id,
      raw_consensus_id = paste0(clonotype_id, "_consensus_2"), stringsAsFactors = FALSE
    )
  }
}

contigs_after_qc <- do.call(rbind, contig_rows)
write.table(contigs_after_qc, file.path(out_dir, "contigs_after_qc.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("wrote contigs_after_qc.tsv: %d rows across %d samples\n",
            nrow(contigs_after_qc), length(samples)))

# ---- synthetic GEX counts: n_genes x n_cells, with a STRONG per-state mean
# shift on distinct marker-gene sets, then a real Seurat pipeline
# (Normalize/ScaleData/PCA/UMAP). CoNGA's own downstream DEG/dotplot step
# (find_gex_cluster_degs -> sc.pl.dotplot) needs louvain GEX clustering to
# find real, separable clusters with real differentially-expressed genes -
# a small number of states (independent of TCR clone count) with big,
# non-overlapping marker sets gives it that; too many/weak states (the
# earlier version scaled state count with n_clones, diluting the signal)
# left too few genes anywhere near significant, and CoNGA's own dotplot call
# crashed on the resulting near-empty gene list (matplotlib
# "left cannot be >= right" from a degenerate GridSpec).
n_genes <- 500
gene_names <- sprintf("Gene%04d", seq_len(n_genes))
n_cells <- length(all_barcodes)
n_gex_states <- 3
markers_per_state <- 100  # broad enough that >=50 genes clear scanpy's HVG
                          # selection - CoNGA's internal PCA hard-requires
                          # at least 50 variable genes (n_components=50);
                          # too few genes at a very strong fold-change
                          # over-concentrates variance into too narrow a set.

base_rate <- rgamma(n_genes, shape = 2, rate = 1)
state_vec <- sample(seq_len(n_gex_states), n_cells, replace = TRUE)
state_markers <- split(sample(seq_len(n_genes), n_gex_states * markers_per_state),
                       rep(seq_len(n_gex_states), each = markers_per_state))

counts <- matrix(0L, nrow = n_genes, ncol = n_cells, dimnames = list(gene_names, all_barcodes))
for (cell in seq_len(n_cells)) {
  lambda <- base_rate
  lambda[state_markers[[state_vec[cell]]]] <- lambda[state_markers[[state_vec[cell]]]] * 6
  counts[, cell] <- rpois(n_genes, lambda * 5)
}
counts <- as(counts, "CsparseMatrix")

rownames(all_meta) <- all_meta$barcode
seu <- CreateSeuratObject(counts = counts, meta.data = all_meta[all_barcodes, , drop = FALSE])

seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = min(200, n_genes), verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 20, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:20, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:20, verbose = FALSE)

saveRDS(seu, file.path(out_dir, "seurat_annotated.rds"))
cat(sprintf("wrote seurat_annotated.rds: %d genes x %d cells, samples: %s\n",
            nrow(seu), ncol(seu), paste(samples, collapse = ", ")))
