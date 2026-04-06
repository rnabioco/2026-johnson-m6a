#!/usr/bin/env Rscript
# Preprocess pileup data for Quarto website
#
# Reads large pileup bed.gz files and produces small summary tables.
# Output: results/hotair/summary/website/

suppressPackageStartupMessages({
  library(data.table)
  library(readr)
  library(dplyr)
})

# --- Paths ---
file_arg <- sub("--file=", "",
  commandArgs(trailingOnly = FALSE)[
    grep("--file=", commandArgs(trailingOnly = FALSE))
  ]
)
if (length(file_arg) > 0) {
  project_dir <- normalizePath(file.path(dirname(file_arg), ".."))
} else {
  project_dir <- normalizePath(".")
}

results_dir <- file.path(project_dir, "results", "hotair")
output_dir  <- file.path(results_dir, "summary", "website")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c("HOTAIR_WT", "HOTAIR_A783U", "HOTAIR_bc06", "HOTAIR_AL")

# Pileup bedMethyl columns (modkit output)
pileup_cols <- c(
  "chrom", "start", "end", "mod_code", "score", "strand",
  "thick_start", "thick_end", "color", "n_valid_cov",
  "pct_modified", "n_mod", "n_canonical", "n_other_mod",
  "n_delete", "n_fail", "n_diff", "n_nocall"
)

# Modification code labels
mod_labels <- c(
  a = "m6A", m = "m5C",
  "17596" = "inosine", "17802" = "pseudouridine",
  "19227" = "Cm", "19228" = "Am", "19229" = "Um", "69426" = "Gm"
)

barcode_to_sample <- c(
  barcode04 = "HOTAIR_WT", barcode05 = "HOTAIR_A783U",
  barcode06 = "HOTAIR_bc06", barcode07 = "HOTAIR_AL"
)

# Helper: parse GENCODE chrom field into metadata columns
# Format: ENST...|ENSG...|OTTHUM...|OTTHUM...|TxName|GeneName|Length|Type|
parse_chrom <- function(dt) {
  dt[, c("transcript_id", "gene_id", "havana_gene", "havana_transcript",
         "transcript_name", "gene_name", "transcript_length", "transcript_type") :=
       tstrsplit(chrom, "\\|", keep = 1:8)]
  dt[, transcript_length := as.integer(transcript_length)]
  dt
}

# Helper: add mod_label column
add_mod_labels <- function(dt) {
  dt[, mod_label := fifelse(mod_code %in% names(mod_labels),
                            mod_labels[mod_code], mod_code)]
  dt
}


# ==========================================================
# 1. Demux summary
# ==========================================================
cat("1. Demux summary\n")
demux_files <- Sys.glob(file.path(
  results_dir, "demux", "read_ids", "*", "demux_summary.tsv.gz"
))
if (length(demux_files) > 0) {
  demux <- read_tsv(demux_files[1], show_col_types = FALSE) |>
    mutate(sample = barcode_to_sample[predicted_barcode]) |>
    arrange(desc(n_reads))
  write_tsv(demux, file.path(output_dir, "demux_summary.tsv"))
  cat("   -> demux_summary.tsv\n")
} else {
  cat("   -> No demux summary found, skipping\n")
}


# ==========================================================
# 2. Alignment stats (samtools flagstat)
# ==========================================================
cat("2. Alignment stats\n")
aln_stats <- lapply(samples, function(s) {
  bam <- file.path(results_dir, "bam", "final", s, paste0(s, ".bam"))
  if (!file.exists(bam)) return(NULL)
  fs <- system(paste("samtools flagstat", bam), intern = TRUE)
  total <- as.numeric(sub(" .*", "", fs[1]))
  mapped_line <- fs[grep("mapped \\(", fs)[1]]
  mapped <- as.numeric(sub(" .*", "", mapped_line))
  mapped_pct <- as.numeric(sub(".*\\(([0-9.]+)%.*", "\\1", mapped_line))
  primary_line <- fs[grep("primary mapped", fs)]
  primary <- if (length(primary_line) > 0) {
    as.numeric(sub(" .*", "", primary_line[1]))
  } else {
    mapped
  }
  data.frame(
    sample = s, total_reads = total, mapped_reads = mapped,
    mapping_rate = mapped_pct, primary_mapped = primary
  )
}) |> bind_rows()

if (nrow(aln_stats) > 0) {
  write_tsv(aln_stats, file.path(output_dir, "alignment_stats.tsv"))
  cat("   -> alignment_stats.tsv\n")
}


# ==========================================================
# 3. HOTAIR pileup (filtered)
# ==========================================================
cat("3. HOTAIR pileup\n")
hotair_list <- lapply(samples, function(s) {
  bed <- file.path(results_dir, "summary", "modkit", s, paste0(s, ".pileup.bed.gz"))
  if (!file.exists(bed)) return(NULL)
  cat("   Reading HOTAIR rows from", s, "\n")
  cmd <- paste0("zcat ", bed, " | grep '|HOTAIR|'")
  dt <- tryCatch(
    fread(cmd = cmd, header = FALSE, sep = "\t", col.names = pileup_cols),
    error = function(e) { cat("   Warning:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(dt) && nrow(dt) > 0) dt[, sample := s]
  dt
})

hotair_all <- rbindlist(Filter(Negate(is.null), hotair_list))
if (nrow(hotair_all) > 0) {
  parse_chrom(hotair_all)
  add_mod_labels(hotair_all)
  fwrite(hotair_all, file.path(output_dir, "hotair_pileup.tsv.gz"),
         sep = "\t", compress = "gzip")
  cat("   -> hotair_pileup.tsv.gz (", nrow(hotair_all), " rows)\n")
}


# ==========================================================
# 4. Top transcripts by coverage
# ==========================================================
cat("4. Top transcripts\n")
top_list <- lapply(samples, function(s) {
  bed <- file.path(results_dir, "summary", "modkit", s, paste0(s, ".pileup.bed.gz"))
  if (!file.exists(bed)) return(NULL)
  cat("   Aggregating", s, "\n")
  cmd <- paste0(
    "zcat ", bed,
    " | awk -F'\\t' '{c=$10+0; if(c>m[$1]) m[$1]=c} END{for(k in m) print k\"\\t\"m[k]}'"
  )
  dt <- tryCatch(
    fread(cmd = cmd, header = FALSE, sep = "\t",
          col.names = c("chrom", "max_coverage")),
    error = function(e) { cat("   Warning:", conditionMessage(e), "\n"); NULL }
  )
  if (!is.null(dt) && nrow(dt) > 0) dt[, sample := s]
  dt
})

top_all <- rbindlist(Filter(Negate(is.null), top_list))
if (nrow(top_all) > 0) {
  top_transcripts <- top_all[, .(
    mean_max_coverage = mean(max_coverage),
    max_coverage = max(max_coverage),
    n_samples = .N
  ), by = chrom][order(-mean_max_coverage)][1:20]

  parse_chrom(top_transcripts)
  write_tsv(as.data.frame(top_transcripts),
            file.path(output_dir, "top_transcripts.tsv"))
  cat("   -> top_transcripts.tsv\n")

  # ==========================================================
  # 5. Mod summary for top transcripts
  # ==========================================================
  cat("5. Mod summary for top transcripts\n")
  top_chroms <- top_transcripts$chrom

  mod_list <- lapply(samples, function(s) {
    bed <- file.path(results_dir, "summary", "modkit", s,
                     paste0(s, ".pileup.bed.gz"))
    if (!file.exists(bed)) return(NULL)
    cat("   Reading top transcripts from", s, "\n")
    # Write patterns to temp file for grep -F
    pf <- tempfile()
    writeLines(top_chroms, pf)
    cmd <- paste0("zcat ", bed, " | grep -F -f ", pf)
    dt <- tryCatch(
      fread(cmd = cmd, header = FALSE, sep = "\t", col.names = pileup_cols),
      error = function(e) { cat("   Warning:", conditionMessage(e), "\n"); NULL }
    )
    unlink(pf)
    if (is.null(dt) || nrow(dt) == 0) return(NULL)
    # Keep only exact matches
    dt <- dt[chrom %in% top_chroms]
    if (nrow(dt) == 0) return(NULL)
    dt[, .(
      mean_pct_modified  = mean(pct_modified, na.rm = TRUE),
      max_pct_modified   = max(pct_modified, na.rm = TRUE),
      mean_coverage      = mean(n_valid_cov, na.rm = TRUE),
      total_positions    = .N,
      positions_with_mod = sum(n_mod > 0)
    ), by = .(chrom, mod_code)][, sample := s]
  })

  mod_summary <- rbindlist(Filter(Negate(is.null), mod_list))
  if (nrow(mod_summary) > 0) {
    parse_chrom(mod_summary)
    add_mod_labels(mod_summary)
    write_tsv(as.data.frame(mod_summary),
              file.path(output_dir, "mod_summary_by_transcript.tsv"))
    cat("   -> mod_summary_by_transcript.tsv\n")
  }
}

cat("\nDone! Output in:", output_dir, "\n")
