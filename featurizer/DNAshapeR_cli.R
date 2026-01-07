#!/usr/bin/env Rscript
# =============================================================================
# DNAshapeR_cli.R
#
# Purpose
# -------
# Compute DNA shape features from an input FASTA using the Bioconductor package
# DNAshapeR and write them to a CSV. This script is designed to be used as a
# backend for the Python featurizer (feature name: "dna_shape"), but can also be
# run standalone from the command line.
#
# Why this exists
# ---------------
# The original featurizer implementation in some pipelines relies on an online
# DNAshape server (RohsLab) to compute shape features. That approach is not
# scalable and may fail due to server availability or rate limits. DNAshapeR
# performs the same style of DNA shape prediction locally and is suitable for
# batch/HPC workflows (including Singularity).
#
# Output modes
# ------------
# This script supports two output modes:
#
#   --mode mean   (recommended for ML featurization)
#       Produces a compact per-sequence table with 4 columns:
#         HelT, MGW, ProT, Roll
#       Each value is the mean across positions within a sequence.
#       Output file defaults to: dna_shape.csv
#
#   --mode full   (high-dimensional; use only if you really want positional features)
#       Produces a wide per-position matrix including multiple shape types
#       (inter-bp and intra-bp shapes). Columns are labeled like:
#         HelT.nt1, HelT.nt2, ..., MGW.nt1, ...
#       This file can be very large.
#
# Inputs / assumptions
# --------------------
# - Input must be a FASTA file with unique sequence headers.
# - Row names in the output CSV correspond to FASTA headers (sequence IDs).
# - DNAshapeR creates intermediate "proc files" next to the FASTA; we remove
#   them by default (unless --keep-proc is provided).
#
# Usage examples
# --------------
# 1) Compact per-sequence output (default):
#    Rscript DNAshapeR_cli.R --fasta locs.fasta --out dna_shape.csv --mode mean
#
# 2) Full per-position output:
#    Rscript DNAshapeR_cli.R --fasta locs.fasta --out DNAshapeR_full.csv --mode full
#
# 3) Keep DNAshapeR intermediate proc files for debugging:
#    Rscript DNAshapeR_cli.R --fasta locs.fasta --out dna_shape.csv --mode mean --keep-proc
#
# Integration with Python featurizer
# ----------------------------------
# Typically your Python wrapper run_dna_shape() will:
#   1) copy/symlink working_dir/locs.fasta into a subfolder
#   2) call this script with --mode mean
#   3) reindex the CSV to match featurizer's canonical loc ordering
#
# Exit status
# -----------
# - Non-zero exit if input FASTA is missing or arguments are invalid.
# =============================================================================

suppressPackageStartupMessages({
  library(DNAshapeR)   # Local DNA shape predictions (MGW, HelT, ProT, Roll, etc.)
  library(Biostrings)  # FASTA IO: readDNAStringSet()
})

# -----------------------------------------------------------------------------
# CLI argument parsing (no external packages required)
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

# Get the value after a flag, e.g. --fasta FILE
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop(paste("Missing value after", flag))
  return(args[idx + 1])
}

# Check whether a flag is present, e.g. --keep-proc
has_flag <- function(flag) any(args == flag)

# Inputs
in_fasta <- get_arg("--fasta", default = "locs.fasta")
out_csv  <- get_arg("--out",   default = "dna_shape.csv")

# Output mode: "mean" (compact per-sequence) or "full" (per-position wide matrix)
mode <- get_arg("--mode", default = "full")  # mean|full

if (!file.exists(in_fasta)) {
  stop(paste("Input FASTA not found:", in_fasta))
}

# Ensure output directory exists
dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)

# Read FASTA to capture sequence names (headers), which we use as rownames
fa <- readDNAStringSet(in_fasta)
seq_names <- names(fa)

# Helper: replace NA/NaN with 0 in vectors (rare, but can happen if all positions NA)
nan0 <- function(x) {
  x[is.na(x) | is.nan(x)] <- 0
  x
}

# -----------------------------------------------------------------------------
# Main computation
# -----------------------------------------------------------------------------
if (mode == "mean") {
  # "mean" mode: match common ML featurization practice: summarize each shape
  # by a single number per sequence (mean across positions).

  # These 4 are the classic minimal DNAshape set used widely:
  # - HelT: Helix twist
  # - MGW: Minor groove width
  # - ProT: Propeller twist
  # - Roll: Roll angle
  shape_types <- c("HelT", "MGW", "ProT", "Roll")

  # getShape() accepts the FASTA path and returns a list of matrices
  # Each matrix is: rows = sequences, cols = positions
  pred <- getShape(in_fasta, shapeType = shape_types)

  # Compute per-sequence mean for each shape type
  means <- lapply(shape_types, function(t) {
    X <- pred[[t]]
    rownames(X) <- seq_names
    m <- rowMeans(X, na.rm = TRUE)
    nan0(m)
  })

  df <- data.frame(means)
  colnames(df) <- shape_types
  rownames(df) <- seq_names

  write.csv(df, out_csv, quote = FALSE)

} else if (mode == "full") {
  # "full" mode: produce a wide, per-position feature matrix.
  # This is much higher dimensional and can be very large.

  # Inter-bp and intra-bp shape types supported by DNAshapeR
  shape_inter <- c("HelT", "Rise", "Roll", "Shift", "Slide", "Tilt")
  shape_intra <- c("Buckle", "Opening", "ProT", "Shear", "Stagger", "Stretch", "MGW")

  pred_inter <- getShape(in_fasta, shapeType = shape_inter)
  pred_intra <- getShape(in_fasta, shapeType = shape_intra)

  # Convert each shape matrix into a data.frame-like matrix with:
  # - NA replaced by 0
  # - rownames as sequence IDs
  # - informative per-position column names (Shape.nt1 ... Shape.ntN)
  fix_mat <- function(X, nm) {
    X[is.na(X)] <- 0
    rownames(X) <- seq_names
    colnames(X) <- paste0(nm, ".nt", seq_len(ncol(X)))
    X
  }

  mats <- c(
    lapply(names(pred_inter), function(nm) fix_mat(pred_inter[[nm]], nm)),
    lapply(names(pred_intra), function(nm) fix_mat(pred_intra[[nm]], nm))
  )

  full <- Reduce(cbind, mats)
  write.csv(full, out_csv, quote = FALSE)

} else {
  stop("Invalid --mode. Use: mean|full")
}

# -----------------------------------------------------------------------------
# Cleanup: remove DNAshapeR intermediate files produced next to the FASTA
# -----------------------------------------------------------------------------
# DNAshapeR generates intermediate "proc files" with names like:
#   <input_fasta_basename>..HelT, <input_fasta_basename>..MGW, etc.
#
# In batch pipelines, leaving these can clutter working directories and consume
# storage. We remove them unless --keep-proc is provided.
base <- basename(in_fasta)

# Escape '.' for regex matching so locs.fasta matches literally
base_escaped <- gsub("([.])", "\\\\1", base)

proc_files <- list.files(
  dirname(in_fasta),
  pattern = paste0("^", base_escaped, "\\.\\..*"),
  full.names = TRUE
)

if (!has_flag("--keep-proc") && length(proc_files) > 0) {
  file.remove(proc_files)
}
