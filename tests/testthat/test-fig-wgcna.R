sim_norm <- function(n_samp = 30, seed = 11) {
  set.seed(seed)
  f1 <- rnorm(n_samp); f2 <- rnorm(n_samp)
  block <- function(f, k) t(sapply(seq_len(k),
    function(i) f * runif(1, 0.8, 1.2) + rnorm(length(f), 0, 0.25)))
  m <- rbind(block(f1, 40), block(f2, 40), matrix(rnorm(70 * n_samp), nrow = 70))
  rownames(m) <- paste0("g", seq_len(nrow(m)))
  colnames(m) <- paste0("s", seq_len(n_samp))
  m
}

build_wg <- function() {
  d <- wgcna_datexpr(sim_norm(), n_genes = 150)
  suppressWarnings(run_wgcna(d, power = 6, min_module_size = 20))
}

test_that("fig_soft_threshold returns a ggplot", {
  skip_if_not_installed("WGCNA")
  sft <- suppressWarnings(wgcna_pick_power(wgcna_datexpr(sim_norm(), 150),
                                           powers = 1:12))
  expect_s3_class(fig_soft_threshold(sft), "ggplot")
})

test_that("fig_module_trait returns a ggplot", {
  skip_if_not_installed("WGCNA")
  wg <- build_wg()
  meta <- data.frame(sample = rownames(wg$MEs),
                     grp = rep(c("A", "B"), length.out = nrow(wg$MEs)))
  tr <- build_traits(meta, rownames(wg$MEs))
  mt <- module_trait_cor(wg$MEs, tr)
  expect_s3_class(fig_module_trait(mt), "ggplot")
  expect_s3_class(fig_module_trait(mt, mode = "publication"), "ggplot")
})

test_that("fig_module_sizes returns a ggplot", {
  skip_if_not_installed("WGCNA")
  wg <- build_wg()
  expect_s3_class(fig_module_sizes(wg), "ggplot")
  expect_s3_class(fig_module_sizes(wg, include_grey = TRUE), "ggplot")
})

test_that("fig_eigengene returns a ggplot with and without groups", {
  skip_if_not_installed("WGCNA")
  wg <- build_wg()
  mod <- setdiff(unique(wg$modules), "grey")[1]
  expect_s3_class(fig_eigengene(wg, mod), "ggplot")
  groups <- stats::setNames(rep(c("A", "B"), length.out = nrow(wg$MEs)),
                            rownames(wg$MEs))
  expect_s3_class(fig_eigengene(wg, mod, groups = groups), "ggplot")
  expect_error(fig_eigengene(wg, "notacolor"), "no eigengene")
})
