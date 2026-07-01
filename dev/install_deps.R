# One-shot dependency installer for RNAflow (cross-platform: Windows / macOS / Linux)
# --------------------------------------------------------------------------------
# Reads DESCRIPTION and installs every Depends / Imports / Suggests package
# (CRAN *and* Bioconductor) via BiocManager, which resolves both repositories.
# Run this once per machine when setting up the dev environment.
#
# Usage (from the package root, i.e. the folder containing DESCRIPTION):
#   source("dev/install_deps.R")
#
# On Windows you also need Rtools (matching your R version, e.g. Rtools44 for
# R 4.4.x) so packages that need compiling can build. Install it from
# https://cran.r-project.org/bin/windows/Rtools/ before running this.

if (!file.exists("DESCRIPTION")) {
  stop("Run this from the RNAflow package root (the folder with DESCRIPTION). ",
       "In R:  setwd('path/to/RNAflow'); source('dev/install_deps.R')",
       call. = FALSE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

# Collect package names from Depends + Imports + Suggests.
desc   <- read.dcf("DESCRIPTION")
fields <- intersect(c("Depends", "Imports", "Suggests"), colnames(desc))
raw    <- paste(desc[, fields], collapse = ",")
pkgs   <- trimws(strsplit(raw, ",")[[1]])
pkgs   <- sub("\\s*\\(.*\\)", "", pkgs)          # strip "(>= x.y.z)" constraints
pkgs   <- unique(pkgs[nzchar(pkgs) & pkgs != "R"])  # drop the "R (>= ...)" entry

installed <- vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)
missing   <- pkgs[!installed]

message(sprintf("RNAflow has %d dependencies; %d already installed, %d missing.",
                length(pkgs), sum(installed), length(missing)))

if (length(missing) > 0) {
  message("Installing: ", paste(missing, collapse = ", "))
  BiocManager::install(missing, update = FALSE, ask = FALSE)
} else {
  message("Everything is already installed.")
}

# devtools drives the dev workflow (load_all / test / document / check).
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}

message("\nDone. Next:\n",
        "  devtools::load_all()   # load the package\n",
        "  devtools::test()       # run the test suite\n",
        "  RNAflow::run_app()     # launch the app")
