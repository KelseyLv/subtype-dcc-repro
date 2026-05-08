# NEMO (Shamir-Lab GitHub): list of matrices features x samples.
# Usage:
#   Rscript run_nemo.R KIRC "D:/.../fea/KIRC" "D:/.../repro_out/KIRC.nemo.tsv"
#
# num.neighbors defaults to max(5, floor(n_samples / num.clusters)) per NEMO paper hint.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_nemo.R <CANCER> <fea_dir> <out.tsv>")
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

suppressPackageStartupMessages(library(NEMO))

omics_fs <- load_omics_fs(fea_dir)
sample_ids <- colnames(omics_fs$rna)
omics_list <- list(omics_fs$CN, omics_fs$meth, omics_fs$miRNA, omics_fs$rna)

n <- length(sample_ids)
num_nei <- max(5L, as.integer(floor(n / k_clust)))

labels <- nemo.clustering(omics_list, num.clusters = k_clust, num.neighbors = num_nei)
labels <- as.integer(labels)
if (is.null(names(labels))) names(labels) <- sample_ids
# ensure order matches sample_ids
labels <- labels[sample_ids]

write_cluster_tsv(sample_ids, labels, out_tsv)
message("Wrote ", out_tsv, " (n=", length(labels), ", k=", k_clust, ", num.neighbors=", num_nei, ")")
