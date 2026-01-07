#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(tibble)
})

# ============================================================
# Feature summarizer
# - Reads per-module feature tables from a featurizer output dir
# - Standardizes the ID column to `sequence_name`
# - Drops known non-feature / metadata columns (e.g., index, name_meta*)
# - Produces:
#     1) summary_<allele>.csv  (wide merged feature table)
#     2) FeatureCategory_<allele>.csv (feature -> category map)
# ============================================================

# ----------------------------
# CLI args
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx == length(args)) stop(paste("Missing value after", flag))
  return(args[idx + 1])
}

dir_input <- get_arg("--dir", NULL)
tfdb_path <- get_arg("--tfdb", NULL)
allele    <- get_arg("--allele", "UNK")
topq      <- as.numeric(get_arg("--topq", "0.9"))

if (is.null(dir_input) || !dir.exists(dir_input)) {
  stop("Provide --dir <features_directory> (must exist).")
}

out_summary <- get_arg("--out-summary",  file.path(dir_input, sprintf("summary_%s.csv", allele)))
out_cat     <- get_arg("--out-category", file.path(dir_input, sprintf("FeatureCategory_%s.csv", allele)))

# ----------------------------
# Helpers
# ----------------------------
read_table_auto <- function(path) {
  if (grepl("\\.gz$", path)) {
    fread(cmd = paste("zcat", shQuote(path)), data.table = FALSE)
  } else {
    fread(path, data.table = FALSE)
  }
}

ensure_unique_names <- function(df) {
  colnames(df) <- make.unique(colnames(df), sep = ".")
  df
}

# Normalize ID column to "sequence_name"
# Handles:
# - sequence_name present
# - one or multiple "name" columns: first is ID, others -> name_meta*
# - blank / X / V1 first column (rownames-like)
standardize_sequence_name <- function(df) {
  df <- as.data.frame(df)

  if ("sequence_name" %in% colnames(df)) {
    return(ensure_unique_names(df))
  }

  if ("name" %in% colnames(df)) {
    idx <- which(colnames(df) == "name")
    colnames(df)[idx[1]] <- "sequence_name"
    if (length(idx) > 1) {
      meta_names <- paste0("name_meta", seq_len(length(idx) - 1))
      colnames(df)[idx[-1]] <- meta_names
    }
    return(ensure_unique_names(df))
  }

  cn <- colnames(df)
  if (length(cn) >= 1 && (cn[1] %in% c("", "X", "V1") || grepl("^X\\.", cn[1]))) {
    colnames(df)[1] <- "sequence_name"
    return(ensure_unique_names(df))
  }

  colnames(df)[1] <- "sequence_name"
  ensure_unique_names(df)
}

# Drop non-feature / metadata columns so they never enter summary/diff
drop_nonfeature_cols <- function(df) {
  df <- as.data.frame(df)
  df <- standardize_sequence_name(df)
  df <- ensure_unique_names(df)

  drop_fixed <- c("index", "name", "name.1", "X", "V1", "")
  drop_regex <- "^name_meta[0-9]+$"

  drop_cols <- unique(c(
    intersect(drop_fixed, colnames(df)),
    grep(drop_regex, colnames(df), value = TRUE)
  ))
  drop_cols <- setdiff(drop_cols, "sequence_name")

  if (length(drop_cols) > 0) {
    df <- df[, setdiff(colnames(df), drop_cols), drop = FALSE]
  }
  df
}

# Robust numeric matrix conversion:
# - coerces each column to numeric (suppresses warnings)
# - sets NAs to 0
# - preserves rownames
as_numeric_matrix <- function(df, id_col = "sequence_name", verbose = TRUE) {
  df <- standardize_sequence_name(df)
  df <- ensure_unique_names(df)
  if (!id_col %in% colnames(df)) stop(sprintf("ID column '%s' not found.", id_col))

  df2 <- df %>% column_to_rownames(id_col)

  bad_cols <- character(0)
  mat <- sapply(colnames(df2), function(cn) {
    orig <- df2[[cn]]
    v <- suppressWarnings(as.numeric(orig))
    # detect non-numeric coercion (exclude genuine NAs)
    if (any(is.na(v)) && any(!is.na(orig))) {
      bad_cols <<- union(bad_cols, cn)
    }
    v
  })

  mat <- as.matrix(mat)
  rownames(mat) <- rownames(df2)
  mat[is.na(mat)] <- 0

  if (verbose && length(bad_cols) > 0) {
    message(sprintf(
      "Note: non-numeric coercion in %d columns (set to 0 where NA): %s",
      length(bad_cols), paste(head(bad_cols, 20), collapse = ", ")
    ))
    if (length(bad_cols) > 20) message("  (truncated; more columns affected)")
  }
  mat
}

safe_full_join <- function(dfs) Reduce(function(x, y) full_join(x, y, by = "sequence_name"), dfs)

fill_numeric_na0 <- function(df) {
  df %>% mutate(across(where(is.numeric), ~ tidyr::replace_na(.x, 0)))
}

# Build feature -> category table
feature_category <- function(dfs, cats,
                             id_cols = c("sequence_name", "sequence_id", "Variant", "variant",
                                         "name", "index", "name.1", "X", "V1",
                                         "name_meta1", "name_meta2", "name_meta3",
                                         "y", "label")) {
  stopifnot(length(dfs) == length(cats))
  map2_dfr(dfs, cats, function(df, cat) {
    cn <- colnames(as.data.frame(df))
    tibble(feature = setdiff(cn, id_cols), category = cat)
  }) %>%
    distinct() %>%
    arrange(category, feature)
}

# ----------------------------
# Expected files in dir_input
# ----------------------------
p_5mer     <- file.path(dir_input, "5mer.csv")
p_deepbind <- file.path(dir_input, "deepbind.csv")
p_poly     <- file.path(dir_input, "polyA_polyT_GC.csv")
p_fimo_sum <- file.path(dir_input, "fimo_summary.csv")

# Optional
p_dnashape    <- file.path(dir_input, "dna_shape.csv")
p_sei         <- file.path(dir_input, "sei.csv")
p_fimo_encode <- file.path(dir_input, "fimo_encode", "fimo.txt")
p_fimo_hg19   <- file.path(dir_input, "fimo_hg19", "fimo.txt")

req <- c(p_5mer, p_deepbind, p_poly, p_fimo_sum)
missing <- req[!file.exists(req)]
if (length(missing) > 0) {
  stop(paste("Missing required files:\n", paste(missing, collapse = "\n")))
}

TFDB <- NULL
if (!is.null(tfdb_path) && file.exists(tfdb_path)) {
  TFDB <- read.delim(tfdb_path) %>% select(Symbol, Family, Family.main)
} else {
  message("TFDB not provided/found; TF-family aggregation will be skipped.")
}

# ----------------------------
# Read base feature tables (standardize + drop metadata)
# ----------------------------
kmer5        <- read_table_auto(p_5mer)     %>% standardize_sequence_name() %>% drop_nonfeature_cols()
deepbind     <- read_table_auto(p_deepbind) %>% standardize_sequence_name() %>% drop_nonfeature_cols()
polyApolyT   <- read_table_auto(p_poly)     %>% standardize_sequence_name() %>% drop_nonfeature_cols()
fimo_summary <- read_table_auto(p_fimo_sum) %>% standardize_sequence_name() %>% drop_nonfeature_cols()

dna_shape <- NULL
if (file.exists(p_dnashape)) {
  dna_shape <- read_table_auto(p_dnashape) %>% standardize_sequence_name() %>% drop_nonfeature_cols()
}

sei <- NULL
if (file.exists(p_sei)) {
  sei <- read_table_auto(p_sei) %>% standardize_sequence_name() %>% drop_nonfeature_cols()
}

# ----------------------------
# Derived features
# ----------------------------
deepbind_mat <- as_numeric_matrix(deepbind, verbose = TRUE)
cutoff_deepbind <- as.numeric(quantile(deepbind_mat, probs = topq, na.rm = TRUE))
deepbind_top <- rowSums(deepbind_mat >= cutoff_deepbind)
deepbind_top <- data.frame(sequence_name = names(deepbind_top), n.deepbind_top = deepbind_top)

kmer5_mat <- as_numeric_matrix(kmer5, verbose = TRUE)
n5 <- rowSums(kmer5_mat != 0)
n5 <- data.frame(sequence_name = names(n5), n.5mer = n5)

poly_small <- polyApolyT %>%
  {
    cn <- colnames(.)
    if (all(c("polyA", "polyT") %in% cn)) rename(., n.polyA = polyA, n.polyT = polyT) else .
  } %>%
  select(-any_of("GC"))

sei_top <- NULL
if (!is.null(sei)) {
  sei_mat <- as_numeric_matrix(sei, verbose = TRUE)
  cutoff_sei <- as.numeric(quantile(sei_mat, probs = topq, na.rm = TRUE))
  sei_top_vec <- rowSums(sei_mat >= cutoff_sei)
  sei_top <- data.frame(sequence_name = names(sei_top_vec), n.sei_top = sei_top_vec)
}

# ----------------------------
# Optional FIMO expansions
# ----------------------------
fimo_encode_hits <- NULL
fimo_hg19_hits <- NULL
fam_encode <- NULL
fam_hg19 <- NULL

if (file.exists(p_fimo_encode)) {
  fimo_encode <- read.delim(p_fimo_encode)
  colnames(fimo_encode)[3] <- "sequence_name"
  if ("motif_alt_id" %in% colnames(fimo_encode)) fimo_encode <- fimo_encode %>% select(-motif_alt_id)
  fimo_encode <- fimo_encode[complete.cases(fimo_encode), , drop = FALSE]
  fimo_encode$Symbol <- gsub("\\_.*", "", fimo_encode$motif_id)

  if (!is.null(TFDB)) {
    fimo_encode <- merge(fimo_encode, TFDB, by = "Symbol", all.x = TRUE)
    fimo_encode$Family.main[is.na(fimo_encode$Family.main)] <- "unknown"
    fam_encode <- fimo_encode %>%
      select(sequence_name, Family.main) %>%
      table() %>%
      as.data.frame.matrix()
    colnames(fam_encode) <- paste0("n.", colnames(fam_encode), "_encode")
    fam_encode <- data.frame(sequence_name = rownames(fam_encode), fam_encode)
  }

  fimo_encode_hits <- fimo_encode %>%
    count(sequence_name, Symbol) %>%
    pivot_wider(
      id_cols = sequence_name,
      names_from = Symbol,
      values_from = n,
      values_fill = 0,
      names_prefix = "fimo_encode_target_gene_hits_"
    )
}

if (file.exists(p_fimo_hg19)) {
  fimo_hg19 <- read.delim(p_fimo_hg19)
  colnames(fimo_hg19)[3] <- "sequence_name"
  if ("motif_alt_id" %in% colnames(fimo_hg19)) fimo_hg19 <- fimo_hg19 %>% select(-motif_alt_id)
  fimo_hg19 <- fimo_hg19[complete.cases(fimo_hg19), , drop = FALSE]
  fimo_hg19$Symbol <- gsub("\\_.*", "", fimo_hg19$motif_id)

  if (!is.null(TFDB)) {
    fimo_hg19 <- merge(fimo_hg19, TFDB, by = "Symbol", all.x = TRUE)
    fimo_hg19$Family.main[is.na(fimo_hg19$Family.main)] <- "unknown"
    fam_hg19 <- fimo_hg19 %>%
      select(sequence_name, Family.main) %>%
      table() %>%
      as.data.frame.matrix()
    colnames(fam_hg19) <- paste0("n.", colnames(fam_hg19), "_hg19")
    fam_hg19 <- data.frame(sequence_name = rownames(fam_hg19), fam_hg19)
  }

  fimo_hg19_hits <- fimo_hg19 %>%
    count(sequence_name, Symbol) %>%
    pivot_wider(
      id_cols = sequence_name,
      names_from = Symbol,
      values_from = n,
      values_fill = 0,
      names_prefix = "fimo_hg19_target_gene_hits_"
    )
}

# ----------------------------
# Merge
# ----------------------------
blocks <- list(
  poly_small,
  n5, kmer5,
  dna_shape,
  fimo_summary,
  fam_encode, fam_hg19,
  fimo_encode_hits, fimo_hg19_hits,
  deepbind_top, deepbind,
  sei_top, sei
)
blocks <- blocks[!vapply(blocks, is.null, logical(1))]

summary_df <- safe_full_join(blocks) %>% fill_numeric_na0()

# ----------------------------
# Feature -> category
# ----------------------------
dfs_for_cat <- list(
  poly_small,
  n5, kmer5,
  dna_shape,
  fimo_summary,
  fam_encode, fam_hg19,
  fimo_encode_hits, fimo_hg19_hits,
  deepbind_top, deepbind,
  sei_top, sei
)

cats_for_cat <- c(
  "polyA_polyT_GC",
  "5mer", "5mer",
  "dna_shape",
  "fimo_summary",
  "fimo_family", "fimo_family",
  "fimo_hits", "fimo_hits",
  "deepbind", "deepbind",
  "sei", "sei"
)

keep_idx <- !vapply(dfs_for_cat, is.null, logical(1))
cat_df <- feature_category(dfs_for_cat[keep_idx], cats_for_cat[keep_idx])

write.csv(summary_df, out_summary, row.names = FALSE)
write.csv(cat_df, out_cat, row.names = FALSE)

message(sprintf("Wrote summary: %s", out_summary))
message(sprintf("Wrote feature categories: %s", out_cat))
message(sprintf("n_rows=%d  n_features=%d", nrow(summary_df), ncol(summary_df) - 1))
