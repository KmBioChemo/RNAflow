# ====================================================================
# RNAflow — developer workflow
# Run these commands from RStudio while developing.
# ====================================================================

# 1. Load the package as if it were installed
devtools::load_all()

# 2. Launch the app for interactive testing
RNAflow::run_app()

# 3. Re-generate documentation (after changing roxygen comments)
devtools::document()

# 4. Run the test suite
devtools::test()

# 5. Full check (CRAN-like)
devtools::check()

# 6. Build the pkgdown site (after first install of pkgdown)
# pkgdown::build_site()

# 7. Generate a demo dataset for quick testing
make_demo <- function(n_genes = 2000, n_samples = 12, seed = 42) {
  set.seed(seed)
  half <- n_samples / 2
  base_mu <- 2^runif(n_genes, 2, 9)
  fc <- rep(1, n_genes)
  fc[1:200] <- 2^runif(200, 0.5, 3)    # 200 up-regulated
  fc[201:400] <- 2^-runif(200, 0.5, 3) # 200 down-regulated
  ctrl <- sapply(seq_len(half), function(i) rpois(n_genes, base_mu))
  trt  <- sapply(seq_len(half), function(i) rpois(n_genes, base_mu * fc))
  counts <- cbind(ctrl, trt)
  rownames(counts) <- paste0("Gene", seq_len(n_genes))
  colnames(counts) <- c(paste0("Ctrl_", seq_len(half)),
                        paste0("Trt_",  seq_len(half)))
  meta <- data.frame(
    sample    = colnames(counts),
    condition = rep(c("Control", "Treatment"), each = half),
    batch     = rep(c("A", "B"), times = half),
    stringsAsFactors = FALSE
  )
  list(counts = counts, metadata = meta)
}

# Write demo CSVs to inst/extdata
write_demo <- function() {
  demo <- make_demo()
  utils::write.csv(
    data.frame(gene = rownames(demo$counts), demo$counts, check.names = FALSE),
    "inst/extdata/demo_counts.csv", row.names = FALSE
  )
  utils::write.csv(demo$metadata, "inst/extdata/demo_metadata.csv", row.names = FALSE)
  message("Wrote inst/extdata/demo_counts.csv and demo_metadata.csv")
}
# write_demo()
