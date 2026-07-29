#!/usr/bin/env Rscript
# Bridge 2+: CLUSTER_TO_SC
# Maps bulk GIANA / GLIPH2 / TCRdist3 cluster assignments onto single cells.
# Adds metadata columns to the Seurat object and writes per-cell TSVs for
# the Consensus Clustering module.
#
# Usage (called from within a Nextflow work directory):
#   Rscript enrich_seurat.R \
#       --seurat_rds    <path>  \
#       --export_cells  <path>  \
#       --tcrdist_radius <int>

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(igraph)
    library(Matrix)
})

# ── Argument parsing ──────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(args, flag, default = NULL) {
    idx <- which(args == flag)
    if (length(idx) == 0) return(default)
    args[idx + 1]
}

seurat_rds_path  <- get_arg(args, "--seurat_rds")
export_cells_path <- get_arg(args, "--export_cells")
tcrdist_radius   <- as.numeric(get_arg(args, "--tcrdist_radius", "24"))

if (is.null(seurat_rds_path) || is.null(export_cells_path)) {
    stop("Usage: enrich_seurat.R --seurat_rds <rds> --export_cells <tsv> [--tcrdist_radius <n>]")
}

# ── Helpers ───────────────────────────────────────────────────────────────────
extract_beta_cdr3 <- function(ctaa) {
    sapply(ctaa, function(x) {
        if (is.na(x) || x == "" || x == "None") return(NA_character_)
        parts <- strsplit(x, ";")[[1]]
        if (length(parts) >= 2) toupper(trimws(parts[2])) else NA_character_
    }, USE.NAMES = FALSE)
}

# ── Load inputs ───────────────────────────────────────────────────────────────
message("[Bridge 2+] Loading Seurat object: ", seurat_rds_path)
seurat_obj <- readRDS(seurat_rds_path)

message("[Bridge 2+] Loading export_cells: ", export_cells_path)
cells <- tryCatch(
    read.delim(export_cells_path, sep = "\t", stringsAsFactors = FALSE),
    error = function(e) stop("Cannot read export_cells: ", e$message)
)

# Determine CDR3b column
if ("CTaa" %in% colnames(cells)) {
    cells$CDR3b <- extract_beta_cdr3(cells$CTaa)
} else if ("junction_aa" %in% colnames(cells)) {
    cells$CDR3b <- toupper(cells$junction_aa)
} else {
    stop("[Bridge 2+] export_cells must contain 'CTaa' or 'junction_aa'")
}

# Determine barcode column
barcode_col <- intersect(c("barcode", "cell_id", "Barcode", "cell", "Cell"),
                         colnames(cells))[1]
if (is.na(barcode_col)) {
    barcode_col <- colnames(cells)[1]
    message("[Bridge 2+] Barcode column not found by name; using first column: ", barcode_col)
}
message("[Bridge 2+] Using barcode column: '", barcode_col, "'")

# ── GIANA cluster mapping ─────────────────────────────────────────────────────
giana_files <- list.files(".", pattern = "_giana\\.txt$", full.names = TRUE)
giana_map   <- NULL

if (length(giana_files) > 0) {
    message("[Bridge 2+] Reading ", length(giana_files), " GIANA file(s)")
    parts <- lapply(giana_files, function(f) {
        tryCatch({
            df <- read.delim(f, sep = "\t", comment.char = "#",
                             stringsAsFactors = FALSE, check.names = FALSE)
            if (ncol(df) < 2) return(NULL)
            cdr3_col    <- colnames(df)[1]
            cluster_col <- if ("cluster" %in% colnames(df)) "cluster" else colnames(df)[2]
            patient     <- sub("_giana\\.txt$", "", basename(f))
            df %>%
                select(CDR3b = all_of(cdr3_col), raw = all_of(cluster_col)) %>%
                mutate(CDR3b         = toupper(CDR3b),
                       giana_cluster = paste0(patient, "_G", raw)) %>%
                select(CDR3b, giana_cluster) %>%
                distinct(CDR3b, .keep_all = TRUE)
        }, error = function(e) {
            message("[Bridge 2+] Warning: skipping GIANA file ", f, ": ", e$message)
            NULL
        })
    })
    giana_map <- bind_rows(Filter(Negate(is.null), parts)) %>%
        distinct(CDR3b, .keep_all = TRUE)
    message("[Bridge 2+] GIANA: ", nrow(giana_map), " unique CDR3b sequences mapped")
} else {
    message("[Bridge 2+] No GIANA files found — skipping GIANA annotation")
}

# ── GLIPH2 cluster mapping ────────────────────────────────────────────────────
gliph2_files <- list.files(".", pattern = "cluster_member_details.*\\.txt$",
                            full.names = TRUE, recursive = TRUE)
gliph2_map   <- NULL

if (length(gliph2_files) > 0) {
    message("[Bridge 2+] Reading ", length(gliph2_files), " GLIPH2 file(s)")
    parts <- lapply(gliph2_files, function(f) {
        tryCatch({
            df <- read.delim(f, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
            if (!"CDR3b" %in% colnames(df)) return(NULL)
            # 'tag' is the GLIPH2 cluster/motif identifier
            cluster_col <- intersect(c("tag", "seq_ID", "ultCDR3b"), colnames(df))[1]
            if (is.na(cluster_col)) cluster_col <- colnames(df)[6]
            # Patient identity comes from the filename (${patient}_cluster_member_details.txt);
            # fall back to the parent directory for legacy layouts.
            patient_dir <- sub("_cluster_member_details.*$", "", basename(f))
            if (patient_dir == basename(f)) patient_dir <- basename(dirname(f))
            df %>%
                select(CDR3b = CDR3b, raw = all_of(cluster_col)) %>%
                mutate(CDR3b          = toupper(CDR3b),
                       gliph2_cluster = paste0(patient_dir, "_", raw)) %>%
                select(CDR3b, gliph2_cluster) %>%
                distinct(CDR3b, .keep_all = TRUE)
        }, error = function(e) {
            message("[Bridge 2+] Warning: skipping GLIPH2 file ", f, ": ", e$message)
            NULL
        })
    })
    gliph2_map <- bind_rows(Filter(Negate(is.null), parts)) %>%
        distinct(CDR3b, .keep_all = TRUE)
    message("[Bridge 2+] GLIPH2: ", nrow(gliph2_map), " unique CDR3b sequences mapped")
} else {
    message("[Bridge 2+] No GLIPH2 files found — skipping GLIPH2 annotation")
}

# ── TCRdist3 cluster mapping (connected components at radius) ─────────────────
clone_df_files <- list.files(".", pattern = "_clone_df\\.csv$", full.names = TRUE)
tcrdist_map    <- NULL

if (length(clone_df_files) > 0) {
    message("[Bridge 2+] Processing ", length(clone_df_files),
            " TCRdist3 sample(s) at radius ", tcrdist_radius)

    parts <- lapply(clone_df_files, function(cdf_path) {
        tryCatch({
            sample_name <- sub("_clone_df\\.csv$", "", basename(cdf_path))
            clone_df    <- read.csv(cdf_path, stringsAsFactors = FALSE,
                                    check.names = FALSE)

            # Find matching distance matrix
            mat_hdf5 <- paste0(sample_name, "_distance_matrix.hdf5")
            mat_csv  <- paste0(sample_name, "_distance_matrix.csv")

            if (file.exists(mat_hdf5) && requireNamespace("rhdf5", quietly = TRUE)) {
                data    <- rhdf5::h5read(mat_hdf5, "data")
                indices <- rhdf5::h5read(mat_hdf5, "indices")
                indptr  <- rhdf5::h5read(mat_hdf5, "indptr")
                shape   <- as.integer(rhdf5::h5read(mat_hdf5, "shape"))
                mat <- Matrix::sparseMatrix(
                    i    = as.integer(indices) + 1L,
                    p    = as.integer(indptr),
                    x    = as.numeric(data),
                    dims = shape,
                    repr = "C"
                )
                # Sentinel -1 encodes true zero-distance pairs (stored as -1 to
                # distinguish from structural zeros in the sparse format)
                mat@x[mat@x == -1] <- 0
                adj <- (mat > 0) & (mat <= tcrdist_radius)
            } else if (file.exists(mat_csv)) {
                mat <- as.matrix(read.csv(mat_csv, header = FALSE))
                adj <- (mat > 0) & (mat <= tcrdist_radius)
            } else {
                message("[Bridge 2+] No distance matrix for sample '", sample_name,
                        "' — skipping")
                return(NULL)
            }

            # Connected components → cluster IDs
            g     <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected",
                                                         diag = FALSE)
            comps <- igraph::components(g)
            n     <- nrow(clone_df)

            cdr3_col <- intersect(c("junction_aa", "CDR3b", "cdr3_b_aa",
                                    "sequence_id"), colnames(clone_df))[1]
            if (is.na(cdr3_col)) {
                message("[Bridge 2+] Cannot identify CDR3b column in clone_df for '",
                        sample_name, "' — skipping")
                return(NULL)
            }

            data.frame(
                CDR3b          = toupper(clone_df[[cdr3_col]][seq_len(n)]),
                tcrdist_cluster = paste0(sample_name, "_T",
                                         comps$membership[seq_len(n)]),
                stringsAsFactors = FALSE
            ) %>% distinct(CDR3b, .keep_all = TRUE)

        }, error = function(e) {
            message("[Bridge 2+] Warning: TCRdist3 clustering failed for '",
                    cdf_path, "': ", e$message)
            NULL
        })
    })

    valid_parts <- Filter(Negate(is.null), parts)
    if (length(valid_parts) == 0) {
        tcrdist_map <- NULL
        message("[Bridge 2+] TCRdist3: all samples failed or had no distance matrix — skipping")
    } else {
        tcrdist_map <- bind_rows(valid_parts) %>%
            distinct(CDR3b, .keep_all = TRUE)
        message("[Bridge 2+] TCRdist3: ", nrow(tcrdist_map), " unique CDR3b sequences mapped")
    }
} else {
    message("[Bridge 2+] No TCRdist3 clone_df files found — skipping TCRdist3 annotation")
}

# ── Write per-cell export TSVs (for Consensus Clustering) ────────────────────
write_export_tsv <- function(cells, cluster_map, cluster_col, out_file) {
    if (is.null(cluster_map) || nrow(cluster_map) == 0) return(invisible(NULL))
    merged   <- left_join(cells, cluster_map, by = "CDR3b")
    assigned <- sum(!is.na(merged[[cluster_col]]))
    message(sprintf("[Bridge 2+] %s: %d / %d cells assigned a %s",
                    out_file, assigned, nrow(merged), cluster_col))
    write.table(merged, out_file, sep = "\t", row.names = FALSE, quote = FALSE)
}

write_export_tsv(cells, giana_map,   "giana_cluster",    "giana_export_cells.tsv")
write_export_tsv(cells, gliph2_map,  "gliph2_cluster",   "gliph2_export_cells.tsv")
write_export_tsv(cells, tcrdist_map, "tcrdist_cluster",  "tcrdist_export_cells.tsv")

# ── Add cluster columns to Seurat metadata ────────────────────────────────────
add_to_seurat <- function(seurat_obj, cells, cluster_map, col_name, barcode_col) {
    if (is.null(cluster_map) || nrow(cluster_map) == 0) return(seurat_obj)
    joined <- left_join(
        cells[, c(barcode_col, "CDR3b"), drop = FALSE],
        cluster_map,
        by = "CDR3b"
    )
    vec     <- setNames(joined[[col_name]], joined[[barcode_col]])
    common  <- intersect(names(vec), colnames(seurat_obj))
    if (length(common) == 0) {
        message("[Bridge 2+] No barcode overlap for '", col_name,
                "' — check barcode format matches Seurat colnames")
        return(seurat_obj)
    }
    meta_vec <- vec[colnames(seurat_obj)]
    seurat_obj[[col_name]] <- unname(meta_vec)
    n_ann <- sum(!is.na(seurat_obj[[col_name]]))
    message(sprintf("[Bridge 2+] Added '%s' to Seurat: %d / %d cells annotated",
                    col_name, n_ann, ncol(seurat_obj)))
    seurat_obj
}

seurat_obj <- add_to_seurat(seurat_obj, cells, giana_map,   "giana_cluster",   barcode_col)
seurat_obj <- add_to_seurat(seurat_obj, cells, gliph2_map,  "gliph2_cluster",  barcode_col)
seurat_obj <- add_to_seurat(seurat_obj, cells, tcrdist_map, "tcrdist_cluster", barcode_col)

# ── Save enriched Seurat ──────────────────────────────────────────────────────
message("[Bridge 2+] Saving enriched Seurat → enriched_seurat.rds")
saveRDS(seurat_obj, "enriched_seurat.rds")
message("[Bridge 2+] Done.")
