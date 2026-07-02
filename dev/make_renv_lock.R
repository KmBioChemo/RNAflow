# Optional: generate an renv.lock pinning the exact CRAN/Bioconductor package
# versions currently installed. Run this deliberately (it activates renv on the
# project, adding renv/ and an .Rprofile that changes how R loads here).
#
# Reproducibility strategy for RNAflow:
#   * PRIMARY  -> the Dockerfile pins the platform (R 4.5 / Bioconductor 3.22 +
#     system libs). This alone is enough for a Bioconductor-heavy app and is the
#     recommended path for sharing / deployment.
#   * OPTIONAL -> renv.lock (this script) pins individual package versions on top
#     of that platform, for byte-for-byte package reproducibility.
#
# Usage (from the package root):
#   Rscript dev/make_renv_lock.R
# Then commit renv.lock. To restore later: renv::restore().
#
# To back out renv afterwards: renv::deactivate() and delete renv/ + .Rprofile.

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

# Discover dependencies from DESCRIPTION + code, snapshot the installed versions.
renv::init(bare = TRUE, restart = FALSE)
renv::snapshot(type = "all", prompt = FALSE)

cat("\nWrote renv.lock. Commit it; restore elsewhere with renv::restore().\n")
