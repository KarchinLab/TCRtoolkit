#!/usr/bin/env Rscript
# Bridge 2+: CLUSTER_TO_SC
# Maps bulk GIANA / GLIPH2 / TCRdist3 cluster assignments onto single cells.
# Adds metadata columns to the Seurat object and writes per-cell TSVs for
# the Consensus Clustering module. The Seurat object is optional: without one
# (bulk data, single-cell VDJ-only route) only the per-cell TSVs are written.
#
# Usage (called from within a Nextflow work directory):
#   Rscript enrich_seurat.R \
#       --seurat_rds    <path>  \
#       --export_cells  <path>  \
#       --tcrdist_radius <int>

suppressPackageStartupMessages({
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

if (is.null(export_cells_path)) {
    stop("Usage: enrich_seurat.R --export_cells <tsv> [--seurat_rds <rds>] [--tcrdist_radius <n>]")
}
# A staged NO_FILE placeholder (0 bytes) means "no GEX object".
has_seurat <- !is.null(seurat_rds_path) && file.exists(seurat_rds_path) &&
    file.size(seurat_rds_path) > 0

# ── Helpers ───────────────────────────────────────────────────────────────────
# scRepertoire's CTaa format is "A:<alpha_cdr3>|B:<beta_cdr3>" (pipe-delimited,
# explicit chain labels) - not the semicolon-delimited, label-free, always-
# alpha-then-beta format this used to assume. Match the "B:" segment by label
# rather than position, and strip the label, so it still works if a cell only
# has one chain or the chains are ordered differently.
extract_beta_cdr3 <- function(ctaa) {
    sapply(ctaa, function(x) {
        if (is.na(x) || x == "" || x == "None") return(NA_character_)
        parts <- strsplit(x, "[|;]")[[1]]
        beta_part <- parts[grepl("^B:", parts)]
        if (length(beta_part) == 0) return(NA_character_)
        toupper(trimws(sub("^B:", "", beta_part[1])))
    }, USE.NAMES = FALSE)
}

# ── Load inputs ───────────────────────────────────────────────────────────────
if (has_seurat) {
    suppressPackageStartupMessages(library(Seurat))
    message("[Bridge 2+] Loading Seurat object: ", seurat_rds_path)
    seurat_obj <- readRDS(seurat_rds_path)
    # Object size, for ranking memory across processes (bin/trace_memory_summary.py).
    cat(sprintf("OBJECT_LOAD\tprocess=CLUSTER_TO_SC\tfile=%s\tfile_mb=%.0f\tobject_mb=%.0f\tcells=%d\tfeatures=%d\tr_max_used_mb=%.0f",
      basename(seurat_rds_path), file.size(seurat_rds_path) / 2^20, as.numeric(utils::object.size(seurat_obj)) / 2^20,
      ncol(seurat_obj), nrow(seurat_obj), sum(gc()[, 6])), "\n", file = stderr())
} else {
    message("[Bridge 2+] No Seurat object - writing per-cell cluster TSVs only")
}

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

# Patient of each row. Every method clusters within a patient, so a cluster labels only
# that patient's cells carrying the CDR3b - not every cell in the cohort that shares it.
patient_col <- intersect(c("patient", "patient_id", "subject_id", "META_PATIENT"), colnames(cells))[1]
if (is.na(patient_col)) {
    message("[Bridge 2+] No patient column in export_cells; matching clusters on CDR3b alone")
    cells$.cluster_patient <- NA_character_
} else {
    message("[Bridge 2+] Matching clusters on (", patient_col, ", CDR3b)")
    cells$.cluster_patient <- as.character(cells[[patient_col]])
}

# Join a (patient, CDR3b, <cluster>) map onto rows of `df`.
join_clusters <- function(df, cluster_map) {
    cluster_map <- dplyr::rename(cluster_map, .cluster_patient = patient)
    if (all(is.na(df$.cluster_patient))) {
        cluster_map <- cluster_map %>% select(-.cluster_patient) %>% distinct(CDR3b, .keep_all = TRUE)
        return(left_join(df, cluster_map, by = "CDR3b"))
    }
    left_join(df, cluster_map, by = c(".cluster_patient", "CDR3b"))
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
                       patient       = patient,
                       giana_cluster = paste0(patient, "_G", raw)) %>%
                select(patient, CDR3b, giana_cluster) %>%
                distinct(CDR3b, .keep_all = TRUE)
        }, error = function(e) {
            message("[Bridge 2+] Warning: skipping GIANA file ", f, ": ", e$message)
            NULL
        })
    })
    giana_map <- bind_rows(Filter(Negate(is.null), parts))
    message("[Bridge 2+] GIANA: ", nrow(giana_map), " (patient, CDR3b) mapped")
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
            # otherwise fall back to the parent directory name.
            patient_dir <- sub("_cluster_member_details.*$", "", basename(f))
            if (patient_dir == basename(f)) patient_dir <- basename(dirname(f))
            df %>%
                select(CDR3b = CDR3b, raw = all_of(cluster_col)) %>%
                mutate(CDR3b          = toupper(CDR3b),
                       patient        = patient_dir,
                       gliph2_cluster = paste0(patient_dir, "_", raw)) %>%
                select(patient, CDR3b, gliph2_cluster) %>%
                distinct(CDR3b, .keep_all = TRUE)
        }, error = function(e) {
            message("[Bridge 2+] Warning: skipping GLIPH2 file ", f, ": ", e$message)
            NULL
        })
    })
    gliph2_map <- bind_rows(Filter(Negate(is.null), parts))
    message("[Bridge 2+] GLIPH2: ", nrow(gliph2_map), " (patient, CDR3b) mapped")
} else {
    message("[Bridge 2+] No GLIPH2 files found — skipping GLIPH2 annotation")
}

# ── TCRdist3 cluster mapping (per patient, precomputed) ───────────────────────
# TCRDIST3_MATRIX cuts each patient's clusters at --tcrdist_radius and writes them as
# <patient>_tcrdist_cluster_map.tsv (one row per clone; tcrdist_cluster empty for clones
# outside any cluster). Clones within the radius are linked, rows sharing a CDR3b are tied
# together, and a connected component with >= 2 distinct CDR3b is a cluster. Reading these
# small tables keeps the raw distance matrices out of this all-patient step.
map_files   <- list.files(".", pattern = "_tcrdist_cluster_map\\.tsv$", full.names = TRUE)
tcrdist_map <- NULL

if (length(map_files) > 0) {
    message("[Bridge 2+] Reading ", length(map_files), " TCRdist3 cluster map(s)")
    parts <- lapply(map_files, function(f) {
        tryCatch({
            m <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE,
                            colClasses = c(tcrdist_cluster = "character"))
            m <- m[!is.na(m$tcrdist_cluster) & nzchar(m$tcrdist_cluster), , drop = FALSE]
            if (!nrow(m)) return(NULL)
            data.frame(patient = as.character(m$patient), CDR3b = toupper(m$cdr3_b_aa),
                       tcrdist_cluster = m$tcrdist_cluster, stringsAsFactors = FALSE) %>%
                distinct(patient, CDR3b, .keep_all = TRUE)
        }, error = function(e) {
            message("[Bridge 2+] Warning: could not read TCRdist3 cluster map '", f, "': ", e$message)
            NULL
        })
    })
    valid_parts <- Filter(Negate(is.null), parts)
    if (length(valid_parts) == 0) {
        message("[Bridge 2+] TCRdist3: no clusters in any patient")
    } else {
        tcrdist_map <- bind_rows(valid_parts)
        message("[Bridge 2+] TCRdist3: ", nrow(tcrdist_map), " clustered (patient, CDR3b) in ",
                dplyr::n_distinct(tcrdist_map$tcrdist_cluster), " clusters (radius ", tcrdist_radius, ")")
    }
} else {
    message("[Bridge 2+] No TCRdist3 cluster maps found — skipping TCRdist3 annotation")
}

# ── Write per-cell export TSVs (for Consensus Clustering) ────────────────────
write_export_tsv <- function(cells, cluster_map, cluster_col, out_file) {
    if (is.null(cluster_map) || nrow(cluster_map) == 0) return(invisible(NULL))
    merged   <- join_clusters(cells, cluster_map) %>% select(-.cluster_patient)
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
    joined <- join_clusters(cells[, c(barcode_col, ".cluster_patient", "CDR3b"), drop = FALSE],
                            cluster_map)
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

if (has_seurat) {
    seurat_obj <- add_to_seurat(seurat_obj, cells, giana_map,   "giana_cluster",   barcode_col)
    seurat_obj <- add_to_seurat(seurat_obj, cells, gliph2_map,  "gliph2_cluster",  barcode_col)
    seurat_obj <- add_to_seurat(seurat_obj, cells, tcrdist_map, "tcrdist_cluster", barcode_col)

    # ── Save enriched Seurat ──────────────────────────────────────────────────
    message("[Bridge 2+] Saving enriched Seurat → enriched_seurat.rds")
    saveRDS(seurat_obj, "enriched_seurat.rds")
}
message("[Bridge 2+] Done.")
