# PINSPlus::SubtypingOmicsData — rows = samples, cols = features per omic.
# Usage:
#   Rscript run_pins.R KIRC "D:/.../fea/KIRC" "D:/.../repro_out/KIRC.pins.tsv"
#
# High-dimensional omics can be slow; ncore>1 speeds PerturbationClustering if available.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript run_pins.R <CANCER> <fea_dir> <out.tsv>")
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

suppressPackageStartupMessages(library(PINSPlus))

omics_fs <- load_omics_fs(fea_dir)
sample_ids <- colnames(omics_fs$rna)

dataList <- list(
  t(omics_fs$CN),
  t(omics_fs$meth),
  t(omics_fs$miRNA),
  t(omics_fs$rna)
)
names(dataList) <- c("CN", "meth", "miRNA", "rna")

ncore <- as.integer(Sys.getenv("PINS_NCORE", unset = "1"))
res <- SubtypingOmicsData(
  dataList = dataList,
  k = k_clust,
  ncore = ncore,
  verbose = FALSE,
  clusteringMethod = "kmeans",
  clusteringOptions = list(nstart = 10)
)

lab <- res$cluster2
labels <- as.integer(factor(lab, levels = sort(unique(as.character(lab)))))

write_cluster_tsv(sample_ids, labels, out_tsv)
message("Wrote ", out_tsv, " (n=", length(labels), ", k=", k_clust, ")")
