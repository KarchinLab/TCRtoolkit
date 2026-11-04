#!/usr/bin/env Rscript
# Pair each VDJ (TCR) sample with its GEX sample, without requiring identical names.
#
# Sourced by VDJ_QC_analysis.qmd and TCell_Integration_Report.qmd (bin/ is on PATH in every
# Nextflow task, locally and on AWS Batch). Base R only, so it runs in any container.
#
# TCR and GEX libraries from the same 10x GEM well are usually named by library type
# ("...-TCR-LIB" vs "...-GEX-LIB"), and finalized names cannot always be changed. Matching is
# therefore tried in order, per VDJ sample:
#   1. exact          - identical name
#   2. normalized     - identical after upper-casing, splitting on - _ . and space, and dropping
#                       library-type tokens (TCR, GEX, LIB, ...) as WHOLE tokens only
#   3. barcode_overlap- no name match: the GEX sample sharing the most cell barcodes, accepted
#                       only if the overlap is at least `min_overlap` AND at least `min_ratio`
#                       times the runner-up (paired libraries share most barcodes; unrelated
#                       ones overlap at chance, typically < 5%)
#   4. unmatched      - none of the above; never guessed
# Barcode overlap is also computed for name matches and reported, so a library mispairing
# (name says A, barcodes say B) is flagged in `note` rather than silently accepted. A name
# match is never overridden by barcodes - the note tells the user to check.
#
# Two VDJ samples may not claim the same GEX sample: the one with the higher barcode overlap
# keeps it, the other becomes unmatched with a note.

DEFAULT_SAMPLE_DROP_TOKENS <- c("TCR", "VDJ", "VDJT", "GEX", "GE", "RNA", "LIB", "LIBRARY")

normalize_sample_name <- function(x, drop_tokens = DEFAULT_SAMPLE_DROP_TOKENS) {
  # accepts c("TCR","GEX") or "TCR,GEX" (the form a Nextflow param arrives in)
  drop <- toupper(trimws(unlist(strsplit(as.character(drop_tokens), ","))))
  vapply(as.character(x), function(s) {
    if (is.na(s)) return(NA_character_)
    tok <- strsplit(toupper(s), "[-_. ]+")[[1]]
    tok <- tok[nzchar(tok) & !(tok %in% drop)]
    if (length(tok) == 0) NA_character_ else paste(tok, collapse = "-")
  }, character(1), USE.NAMES = FALSE)
}

# The 10x cell barcode itself: the last run of >= 12 nucleotides, wherever it sits
# ("ACGT..-1", "SAMPLE_ACGT..-1", "ACGT.._SAMPLE"). Sample prefixes and -N suffixes vary
# between tools, the nucleotide barcode does not.
core_barcode <- function(x) {
  x <- toupper(as.character(x))
  m <- regmatches(x, gregexpr("[ACGTN]{12,}", x))
  vapply(m, function(v) if (length(v)) v[length(v)] else NA_character_, character(1))
}

# vdj: data.frame(sample, barcode)   one row per contig or cell, duplicates fine
# gex: data.frame(sample, barcode)   one row per GEX cell
reconcile_sample_names <- function(vdj, gex,
                                   drop_tokens = DEFAULT_SAMPLE_DROP_TOKENS,
                                   min_overlap = 0.2, min_ratio = 3) {
  stopifnot(all(c("sample", "barcode") %in% names(vdj)),
            all(c("sample", "barcode") %in% names(gex)))
  vdj <- data.frame(sample = as.character(vdj$sample), bc = core_barcode(vdj$barcode),
                    stringsAsFactors = FALSE)
  gex <- data.frame(sample = as.character(gex$sample), bc = core_barcode(gex$barcode),
                    stringsAsFactors = FALSE)
  vdj <- unique(vdj[!is.na(vdj$sample) & !is.na(vdj$bc), ])
  gex <- unique(gex[!is.na(gex$sample) & !is.na(gex$bc), ])

  vdj_samples <- sort(unique(vdj$sample))
  gex_samples <- sort(unique(gex$sample))
  gex_sets <- split(gex$bc, gex$sample)

  # overlap[v, g] = fraction of VDJ sample v's barcodes present in GEX sample g
  ov <- matrix(0, nrow = length(vdj_samples), ncol = length(gex_samples),
               dimnames = list(vdj_samples, gex_samples))
  for (v in vdj_samples) {
    vb <- vdj$bc[vdj$sample == v]
    if (length(vb) == 0) next
    for (g in gex_samples) ov[v, g] <- sum(vb %in% gex_sets[[g]]) / length(vb)
  }

  gex_norm <- normalize_sample_name(gex_samples, drop_tokens)
  norm_count <- table(gex_norm[!is.na(gex_norm)])

  res <- lapply(vdj_samples, function(v) {
    rates <- if (length(gex_samples)) ov[v, ] else numeric(0)
    ord <- order(-rates)
    best_g  <- if (length(ord)) gex_samples[ord[1]] else NA_character_
    best_r  <- if (length(ord)) rates[ord[1]] else 0
    second  <- if (length(ord) > 1) rates[ord[2]] else 0

    gex_sample <- NA_character_; method <- "unmatched"; note <- ""
    vn <- normalize_sample_name(v, drop_tokens)
    if (v %in% gex_samples) {
      gex_sample <- v; method <- "exact"
    } else if (!is.na(vn) && vn %in% names(norm_count)) {
      if (norm_count[[vn]] == 1) {
        gex_sample <- gex_samples[which(gex_norm == vn)]; method <- "normalized"
      } else {
        note <- paste0("normalized name '", vn, "' matches ", norm_count[[vn]],
                       " GEX samples; ambiguous, not used")
      }
    }
    if (is.na(gex_sample) && best_r >= min_overlap &&
        (second == 0 || best_r >= min_ratio * second)) {
      gex_sample <- best_g; method <- "barcode_overlap"
    }
    rate <- if (!is.na(gex_sample)) rates[[gex_sample]] else NA_real_
    if (method %in% c("exact", "normalized") && !is.na(best_g) && best_g != gex_sample &&
        best_r >= min_overlap && best_r >= min_ratio * max(rate, 1e-9)) {
      note <- paste0("possible library mispairing: name matches ", gex_sample,
                     sprintf(" (%.1f%% barcodes)", 100 * rate), " but barcodes match ",
                     best_g, sprintf(" (%.1f%%)", 100 * best_r))
    } else if (method %in% c("exact", "normalized") && rate < min_overlap) {
      note <- sprintf("low barcode overlap with matched GEX sample (%.1f%%)", 100 * rate)
    } else if (method == "unmatched" && !nzchar(note)) {
      note <- if (length(gex_samples) == 0) "no GEX samples" else
        sprintf("no name match; best barcode overlap %s = %.1f%% (needs >= %.0f%% and %gx runner-up)",
                best_g, 100 * best_r, 100 * min_overlap, min_ratio)
    }
    data.frame(vdj_sample = v, gex_sample = gex_sample, method = method,
               vdj_barcodes = sum(vdj$sample == v),
               barcode_overlap = round(rate, 4),
               best_barcode_gex = best_g, best_barcode_overlap = round(best_r, 4),
               note = note, stringsAsFactors = FALSE)
  })
  res <- do.call(rbind, res)
  if (is.null(res)) {
    return(data.frame(vdj_sample = character(), gex_sample = character(), method = character(),
                      vdj_barcodes = integer(), barcode_overlap = numeric(),
                      best_barcode_gex = character(), best_barcode_overlap = numeric(),
                      note = character(), stringsAsFactors = FALSE))
  }

  # one GEX sample per VDJ sample: keep the claimant with the highest barcode overlap
  claimed <- res$gex_sample[!is.na(res$gex_sample)]
  for (g in unique(claimed[duplicated(claimed)])) {
    idx <- which(res$gex_sample == g)
    keep <- idx[which.max(ifelse(is.na(res$barcode_overlap[idx]), -1, res$barcode_overlap[idx]))]
    for (i in setdiff(idx, keep)) {
      res$note[i] <- paste0("also matched ", g, " (", res$method[i], ") but ",
                            res$vdj_sample[keep], " has higher barcode overlap; left unmatched")
      res$gex_sample[i] <- NA_character_; res$method[i] <- "unmatched"
      res$barcode_overlap[i] <- NA_real_
    }
  }
  res
}

# Named vector vdj_sample -> name to use (GEX name when matched, unchanged otherwise).
sample_name_map <- function(recon) {
  setNames(ifelse(is.na(recon$gex_sample), recon$vdj_sample, recon$gex_sample), recon$vdj_sample)
}

print_sample_reconciliation <- function(recon) {
  cat("=== VDJ -> GEX sample pairing ===\n")
  for (i in seq_len(nrow(recon))) {
    r <- recon[i, ]
    cat(sprintf("  %-40s -> %-40s [%s%s]%s\n", r$vdj_sample,
                ifelse(is.na(r$gex_sample), "(none)", r$gex_sample), r$method,
                ifelse(is.na(r$barcode_overlap), "", sprintf(", %.1f%% barcodes", 100 * r$barcode_overlap)),
                ifelse(nzchar(r$note), paste0("\n      NOTE: ", r$note), "")))
  }
  tab <- table(factor(recon$method, levels = c("exact", "normalized", "barcode_overlap", "unmatched")))
  cat(sprintf("  exact=%d normalized=%d barcode_overlap=%d unmatched=%d\n",
              tab[["exact"]], tab[["normalized"]], tab[["barcode_overlap"]], tab[["unmatched"]]))
  invisible(recon)
}
