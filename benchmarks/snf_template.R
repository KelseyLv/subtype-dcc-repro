#!/usr/bin/env Rscript
# =============================================================================
# SNF baseline (CRAN SNFtool) — canonical script for this repo.
# Paper: Wang et al., Nature Methods 2014; Subtype-DCC Table 1 compares SNF.
#
# Prerequisite (once):
#   Rscript benchmarks/install_snf_only.R
#
# Usage (from repo root; paths can be absolute):
#   Rscript benchmarks/snf_template.R KIRC vendor/Subtype-DCC/subtype_file/fea/KIRC repro_out/KIRC.snf.tsv
#
# Output: TSV with columns sample_name, cluster (1..K) for repro/evaluate.py
# Cluster count K matches Subtype-DCC train.py cancer_dict.
# =============================================================================

suppressPackageStartupMessages(library(SNFtool))

BLOCK_FILES <- c(CN = "CN.fea", meth = "meth.fea", miRNA = "miRNA.fea", rna = "rna.fea")
CANCER_K <- c(
  BRCA = 5L, BLCA = 5L, KIRC = 4L, LUAD = 3L, PAAD = 2L,
  SKCM = 4L, STAD = 3L, UCEC = 4L, UVM = 4L
)

read_fea_csv <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  m <- as.matrix(read.csv(path, row.names = 1, check.names = FALSE))
  storage.mode(m) <- "numeric"
  m
}

load_omics_fs <- function(fea_dir) {
  fea_dir <- normalizePath(fea_dir, winslash = "/", mustWork = TRUE)
  mats <- lapply(file.path(fea_dir, BLOCK_FILES), read_fea_csv)
  names(mats) <- names(BLOCK_FILES)
  common <- Reduce(intersect, lapply(mats, colnames))
  if (length(common) < 10) stop("Too few common samples across omics in ", fea_dir)
  for (nm in names(mats)) mats[[nm]] <- mats[[nm]][, common, drop = FALSE]
  mats
}

write_cluster_tsv <- function(sample_ids, labels, out_path) {
  df <- data.frame(sample_name = sample_ids, cluster = as.integer(labels), stringsAsFactors = FALSE)
  out_path <- normalizePath(out_path, winslash = "/", mustWork = FALSE)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, file = out_path, sep = "\t", row.names = FALSE, quote = FALSE)
}

argv <- commandArgs(trailingOnly = TRUE)
if (length(argv) < 3) {
  stop("Usage: Rscript benchmarks/snf_template.R <CANCER> <fea_dir> <out.tsv>\n",
       "Example: Rscript benchmarks/snf_template.R KIRC vendor/Subtype-DCC/subtype_file/fea/KIRC repro_out/KIRC.snf.tsv",
       call. = FALSE)
}

cancer <- toupper(argv[[1]])
fea_dir <- argv[[2]]
out_tsv <- argv[[3]]

if (!cancer %in% names(CANCER_K)) stop("Unknown cancer: ", cancer)
k_clust <- as.integer(CANCER_K[[cancer]])

omics_fs <- load_omics_fs(fea_dir)
sample_ids <- colnames(omics_fs$rna)

# SNF vignette: each view = rows patients, cols features
dataL <- lapply(omics_fs, function(m) standardNormalization(t(m)))

K_nei <- 20L
alpha <- 0.5
T_iter <- 20L

distL <- lapply(dataL, function(x) (dist2(as.matrix(x), as.matrix(x)))^(1 / 2))
affL <- lapply(distL, function(d) affinityMatrix(d, K_nei, alpha))
W <- SNF(affL, K_nei, T_iter)
labels <- as.integer(spectralClustering(W, k_clust))

write_cluster_tsv(sample_ids, labels, out_tsv)
message("Wrote ", out_tsv, " n=", length(labels), " k=", k_clust)
