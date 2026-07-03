make_fake_de <- function(n = 200, seed = 1) {
  set.seed(seed)
  data.frame(
    gene = paste0("g", seq_len(n)),
    baseMean = rgamma(n, 2, 0.1),
    log2FoldChange = rnorm(n, 0, 2),
    lfcSE = runif(n, 0.1, 0.5),
    stat  = rnorm(n),
    pvalue = runif(n),
    padj   = pmin(runif(n) * 0.5, 1),
    stringsAsFactors = FALSE
  )
}

test_that("fig_volcano returns a ggplot", {
  res <- make_fake_de()
  p <- fig_volcano(res, lfc_thr = 1, padj_thr = 0.05, n_label = 5)
  expect_s3_class(p, "ggplot")
})

test_that("fig_volcano works in publication mode", {
  res <- make_fake_de()
  p <- fig_volcano(res, mode = "publication")
  expect_s3_class(p, "ggplot")
})

test_that("fig_volcano_interactive returns a plotly object", {
  res <- make_fake_de()
  p <- fig_volcano_interactive(res, lfc_thr = 1, padj_thr = 0.05)
  expect_s3_class(p, "plotly")
})

test_that("prep_volcano_data adds expected columns", {
  res <- make_fake_de(20)
  df <- prep_volcano_data(res, 1, 0.05)
  expect_true(all(c("lfc", "padj2", "nlog", "reg", "tip") %in% colnames(df)))
  expect_s3_class(df$reg, "factor")
  expect_equal(levels(df$reg), c("Up", "Down", "NS"))
})

test_that("theme_publication / theme_exploration are valid ggplot themes", {
  expect_s3_class(theme_publication(), "theme")
  expect_s3_class(theme_exploration(), "theme")
  expect_s3_class(fig_theme("publication"), "theme")
  expect_s3_class(fig_theme("exploration"), "theme")
})
