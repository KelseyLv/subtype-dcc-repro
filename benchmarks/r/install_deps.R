# Deprecated: one-shot multi-package installs often hit CRAN timeouts on slow networks.
# The previous background run of this script failed for that reason.
#
# Use instead (from repo root):
#   Rscript benchmarks/install_snf_only.R
#   Rscript benchmarks/snf_template.R <CANCER> <fea_dir> <out.tsv>
#
# Optional later: remotes::install_github("Shamir-Lab/NEMO/NEMO"); install.packages("PINSPlus")

message(
  "benchmarks/r/install_deps.R is deprecated.\n",
  "Use benchmarks/install_snf_only.R + benchmarks/snf_template.R (see comments in this file)."
)
invisible(TRUE)
