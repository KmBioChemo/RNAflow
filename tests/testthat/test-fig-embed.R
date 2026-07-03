make_vst_like <- function(n_genes = 200, n_samp = 8, seed = 1) {
  set.seed(seed)
  m <- matrix(rnorm(n_genes * n_samp, 8, 2), nrow = n_genes,
              dimnames = list(paste0("g", seq_len(n_genes)),
                              paste0("s", seq_len(n_samp))))
  # inject group structure so an embedding has something to separate
  m[, 1:(n_samp / 2)] <- m[, 1:(n_samp / 2)] + 3
  m
}

test_that("compute_umap returns a 2D embedding with one row per sample", {
  skip_if_not_installed("uwot")
  m <- make_vst_like()
  out <- compute_umap(m, n_top = 100, n_neighbors = 5)
  expect_named(out, c("scores", "n_used", "n_neighbors"))
  expect_s3_class(out$scores, "data.frame")
  expect_equal(nrow(out$scores), ncol(m))
  expect_setequal(colnames(out$scores), c("UMAP1", "UMAP2", "sample"))
  expect_setequal(out$scores$sample, colnames(m))
  expect_true(all(is.finite(out$scores$UMAP1)))
})

test_that("compute_umap is deterministic and restores the RNG state", {
  skip_if_not_installed("uwot")
  m <- make_vst_like()
  set.seed(999); before <- runif(1)
  set.seed(999)
  a <- compute_umap(m, n_top = 100, n_neighbors = 5, seed = 7)
  after <- runif(1)                     # RNG stream must be unchanged by umap()
  b <- compute_umap(m, n_top = 100, n_neighbors = 5, seed = 7)
  expect_equal(a$scores, b$scores)      # deterministic for a fixed seed
  expect_equal(before, after)           # global RNG restored
})

test_that("compute_umap clamps n_neighbors and rejects tiny sample sets", {
  skip_if_not_installed("uwot")
  m <- make_vst_like(n_samp = 6)
  out <- compute_umap(m, n_top = 50, n_neighbors = 100)   # over-large request
  expect_lte(out$n_neighbors, ncol(m) - 1)
  expect_error(compute_umap(make_vst_like(n_samp = 3)), "at least 4 samples")
})

test_that("fig_umap and fig_pca_3d build plotly objects", {
  skip_if_not_installed("uwot")
  m <- make_vst_like()
  md <- data.frame(sample = colnames(m),
                   group = rep(c("A", "B"), each = ncol(m) / 2))
  expect_s3_class(fig_umap(m, md, n_top = 100, n_neighbors = 5,
                           color_by = "group"), "plotly")
  expect_s3_class(fig_pca_3d(m, md, n_top = 100, color_by = "group"), "plotly")
})

test_that("fig_pca_3d needs at least 4 samples", {
  m <- make_vst_like(n_samp = 3)
  expect_error(fig_pca_3d(m), "at least 4 samples")
})
