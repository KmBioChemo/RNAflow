qc_res <- function(n = 300, seed = 8) {
  set.seed(seed)
  data.frame(
    gene = paste0("g", seq_len(n)), baseMean = rgamma(n, 2, 0.05),
    log2FoldChange = rnorm(n, 0, 2), lfcSE = 0.2, stat = rnorm(n),
    pvalue = runif(n), padj = pmin(runif(n) * 0.4, 1),
    stringsAsFactors = FALSE)
}

qc_norm <- function(g = 200, s = 8) {
  set.seed(3)
  m <- matrix(rnorm(g * s, 8, 2), nrow = g,
              dimnames = list(paste0("g", seq_len(g)), paste0("s", seq_len(s))))
  m
}

test_that("fig_pval_hist returns a ggplot and needs a pvalue column", {
  expect_s3_class(fig_pval_hist(qc_res()), "ggplot")
  expect_error(fig_pval_hist(qc_res()[, c("gene", "padj")]), "pvalue")
})

test_that("fig_ma returns a ggplot", {
  expect_s3_class(fig_ma(qc_res()), "ggplot")
  expect_s3_class(fig_ma(qc_res(), padj_thr = 0.1, mode = "publication"), "ggplot")
})

test_that("fig_sample_cor returns a pheatmap", {
  skip_if_not_installed("pheatmap")
  m <- qc_norm()
  meta <- data.frame(sample = colnames(m), group = rep(c("A", "B"), 4))
  expect_s3_class(fig_sample_cor(m, meta), "pheatmap")
  expect_s3_class(fig_sample_cor(m), "pheatmap")                 # no metadata
  expect_error(fig_sample_cor(m[, 1, drop = FALSE]), "at least 2 samples")
})

test_that("fig_lib_sizes returns a ggplot with and without metadata", {
  counts <- matrix(rpois(200 * 6, 100), nrow = 200,
                   dimnames = list(NULL, paste0("s", 1:6)))
  meta <- data.frame(sample = paste0("s", 1:6), group = rep(c("A", "B"), 3))
  expect_s3_class(fig_lib_sizes(counts, meta), "ggplot")
  expect_s3_class(fig_lib_sizes(counts), "ggplot")
})
