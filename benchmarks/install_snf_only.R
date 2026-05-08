# Minimal CRAN install for SNF baseline (SNFtool only).
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = max(600, getOption("timeout")))
if (!requireNamespace("SNFtool", quietly = TRUE)) {
  install.packages("SNFtool", dependencies = c("Depends", "Imports", "LinkingTo"))
}
message("SNFtool ready.")
