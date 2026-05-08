# Quick load test after install_deps.R
suppressPackageStartupMessages({
  library(SNFtool)
  library(NEMO)
  library(PINSPlus)
})
message("SNFtool, NEMO, PINSPlus: OK")
if (requireNamespace("LRAcluster", quietly = TRUE)) {
  suppressPackageStartupMessages(library(LRAcluster))
  message("LRAcluster: OK")
} else {
  message("LRAcluster: not installed (optional)")
}
message("test_packages.R finished.")
