# Shared FEA readers: Subtype-GAN style CSV, rows = features, cols = samples (barcodes).

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

#' @return named list of matrices (features x samples), same column order as rna.fea
load_omics_fs <- function(fea_dir) {
  fea_dir <- normalizePath(fea_dir, winslash = "/", mustWork = TRUE)
  paths <- file.path(fea_dir, BLOCK_FILES)
  mats <- lapply(paths, read_fea_csv)
  common <- Reduce(intersect, lapply(mats, colnames))
  if (length(common) < 10) {
    stop("Too few common samples across omics in ", fea_dir)
  }
  for (nm in names(mats)) {
    mats[[nm]] <- mats[[nm]][, common, drop = FALSE]
  }
  mats
}

write_cluster_tsv <- function(sample_ids, labels, out_path) {
  if (length(sample_ids) != length(labels)) stop("sample/label length mismatch")
  df <- data.frame(sample_name = sample_ids, cluster = as.integer(labels), stringsAsFactors = FALSE)
  out_path <- normalizePath(out_path, winslash = "/", mustWork = FALSE)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  write.table(df, file = out_path, sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(out_path)
}
