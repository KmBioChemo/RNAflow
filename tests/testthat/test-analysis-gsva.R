make_gsva_data <- function(n_genes = 120, n_samp = 10, seed = 4) {
  set.seed(seed)
  m <- matrix(rnorm(n_genes * n_samp, 8, 2), nrow = n_genes,
              dimnames = list(paste0("g", seq_len(n_genes)),
                              paste0("s", seq_len(n_samp))))
  sets <- lapply(1:6, function(i) paste0("g", sample.int(n_genes, 12)))
  names(sets) <- paste0("SET_", 1:6)
  list(m = m, sets = sets)
}

test_that("run_gsva returns a sets x samples score matrix (gsva + ssgsea)", {
  skip_if_not_installed("GSVA")
  d <- make_gsva_data()
  for (meth in c("gsva", "ssgsea")) {
    es <- run_gsva(d$m, d$sets, method = meth, min_size = 5)
    expect_true(is.matrix(es))
    expect_equal(ncol(es), ncol(d$m))          # one column per sample
    expect_lte(nrow(es), length(d$sets))       # one row per passing set
    expect_setequal(colnames(es), colnames(d$m))
    expect_true(all(is.finite(es)))
  }
})

test_that("run_gsva validates inputs", {
  skip_if_not_installed("GSVA")
  d <- make_gsva_data()
  expect_error(run_gsva(d$m, list()), "non-empty")
  m2 <- d$m; rownames(m2) <- NULL
  expect_error(run_gsva(m2, d$sets), "rownames")
})

test_that("fig_gsva_heatmap builds a pheatmap from scores", {
  skip_if_not_installed("GSVA")
  skip_if_not_installed("pheatmap")
  d <- make_gsva_data()
  es <- run_gsva(d$m, d$sets, method = "gsva", min_size = 5)
  md <- data.frame(sample = colnames(d$m),
                   group = rep(c("A", "B"), length.out = ncol(d$m)))
  p <- fig_gsva_heatmap(es, md, n_top = 6)
  expect_s3_class(p, "pheatmap")
})
