# LRAcluster (Bioconductor): low-rank + clustering on integrated coordinates.
# Usage:
#   Rscript run_lra.R KIRC "D:/.../fea/KIRC" "D:/.../repro_out/KIRC.lra.tsv"
#
# Uses gaussian type for continuous matrices (features x samples).

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_lra.R <CANCER> <fea_dir> <out.tsv>")
}

cancer <- toupper(args[[1]])
fea_dir <- args[[2]]
out_tsv <- args[[3]]

script_dir <- tryCatch({
  ca <- commandArgs(FALSE)
  i <- grep("^--file=", ca, value = TRUE)
  dirname(normalizePath(sub("^--file=", "", i[1])))
}, error = function(e) getwd())

source(file.path(script_dir, "fea_io.R"))

if (!cancer %in% names(CANCER_K)) stop("Unknown cancer: ", cancer)
k_clust <- as.integer(CANCER_K[[cancer]])

suppressPackageStartupMessages(library(LRAcluster))

omics_fs <- load_omics_fs(fea_dir)
sample_ids <- colnames(omics_fs$rna)

data <- list(omics_fs$CN, omics_fs$meth, omics_fs$miRNA, omics_fs$rna)
types <- list("gaussian", "gaussian", "gaussian", "gaussian")

set.seed(21)
fit <- LRAcluster(data = data, types = types, dimension = min(50L, max(3L, ncol(omics_fs$rna) - 1L)))

coord <- fit$coordinate
if (is.null(coord) || nrow(coord) != length(sample_ids)) {
  stop("Unexpected LRAcluster output; check LRAcluster version / inputs.")
}

km <- kmeans(coord, centers = k_clust, nstart = 25, iter.max = 100)
labels <- as.integer(km$cluster)

write_cluster_tsv(sample_ids, labels, out_tsv)
message("Wrote ", out_tsv, " (n=", length(labels), ", k=", k_clust, ")")
