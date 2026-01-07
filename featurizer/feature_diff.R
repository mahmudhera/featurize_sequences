#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default=NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop(paste("Missing value after", flag))
  args[idx + 1]
}

dir_cell <- get_arg("--cell-dir", NULL)   # e.g. .../HEPG2/features
out_path <- get_arg("--out", file.path(dir_cell, "delta_alt_minus_ref.csv"))

if (is.null(dir_cell) || !dir.exists(dir_cell)) stop("Provide --cell-dir")

ref_path <- file.path(dir_cell, "ref", "summary_ref.csv")
alt_path <- file.path(dir_cell, "alt", "summary_alt.csv")

if (!file.exists(ref_path)) stop(paste("Missing:", ref_path))
if (!file.exists(alt_path)) stop(paste("Missing:", alt_path))

ref <- fread(ref_path, data.table=FALSE)
alt <- fread(alt_path, data.table=FALSE)

# Ensure consistent ordering by sequence_name
ref <- ref %>% arrange(sequence_name)
alt <- alt %>% arrange(sequence_name)

if (!all(ref$sequence_name == alt$sequence_name)) {
  stop("sequence_name mismatch between ref and alt summaries")
}

# Compute numeric deltas (alt - ref) for shared columns
common <- intersect(colnames(ref), colnames(alt))
common <- setdiff(common, "sequence_name")

ref_num <- ref[, common, drop=FALSE]
alt_num <- alt[, common, drop=FALSE]

# Convert to numeric safely
for (c in common) {
  ref_num[[c]] <- as.numeric(ref_num[[c]])
  alt_num[[c]] <- as.numeric(alt_num[[c]])
}

delta <- alt_num - ref_num
delta <- cbind(sequence_name = ref$sequence_name, delta)

fwrite(delta, out_path)
cat("Wrote:", out_path, "\n")
