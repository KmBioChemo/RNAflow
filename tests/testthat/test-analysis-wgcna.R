# Small modular dataset: 2 co-expressed blocks driven by latent factors,
# plus noise genes. genes x samples (RNAflow's normalized-matrix convention).
sim_norm <- function(n_samp = 30, seed = 11) {
  set.seed(seed)
  f1 <- rnorm(n_samp); f2 <- rnorm(n_samp)
  block <- function(f, k) t(sapply(seq_len(k),
    function(i) f * runif(1, 0.8, 1.2) + rnorm(length(f), 0, 0.25)))
  m <- rbind(block(f1, 40), block(f2, 40),
             matrix(rnorm(70 * n_samp), nrow = 70))      # noise
  rownames(m) <- paste0("g", seq_len(nrow(m)))
  colnames(m) <- paste0("s", seq_len(n_samp))
  m
}

test_that("wgcna_datexpr selects top-variance genes and transposes", {
  m <- sim_norm()
  d <- wgcna_datexpr(m, n_genes = 100)
  expect_equal(nrow(d), ncol(m))     # samples in rows
  expect_equal(ncol(d), 100)         # genes in cols
  expect_equal(rownames(d), colnames(m))
})

test_that("wgcna_datexpr rejects tiny sample counts", {
  expect_error(wgcna_datexpr(matrix(1:6, nrow = 3)), "at least 4 samples")
})

test_that("build_traits expands annotations into indicator columns", {
  meta <- data.frame(sample = paste0("s", 1:6),
                     cond = rep(c("A", "B", "C"), 2),
                     batch = rep(c("X", "Y"), 3))
  tr <- build_traits(meta, paste0("s", 1:6))
  expect_equal(nrow(tr), 6)
  expect_true(all(c("cond: A", "cond: B", "cond: C",
                    "batch: X", "batch: Y") %in% colnames(tr)))
  expect_true(all(tr %in% c(0, 1)))
})

test_that("wgcna_pick_power returns fit indices and a numeric suggestion", {
  skip_if_not_installed("WGCNA")
  d <- wgcna_datexpr(sim_norm(), n_genes = 150)
  sft <- suppressWarnings(wgcna_pick_power(d, powers = 1:12))
  expect_true(is.data.frame(sft$fit_indices))
  expect_true("SFT.R.sq" %in% colnames(sft$fit_indices))
  expect_true(is.numeric(sft$suggested) && sft$suggested >= 1)
})

test_that("run_wgcna detects modules and downstream helpers work", {
  skip_if_not_installed("WGCNA")
  d <- wgcna_datexpr(sim_norm(), n_genes = 150)
  wg <- suppressWarnings(run_wgcna(d, power = 6, min_module_size = 20))

  expect_named(wg$modules)
  expect_equal(length(wg$modules), ncol(d))
  expect_equal(nrow(wg$MEs), nrow(d))
  # The two planted blocks should yield at least one non-grey module
  expect_gte(length(setdiff(unique(wg$modules), "grey")), 1)

  ms <- module_summary(wg)
  expect_true(all(c("module", "n_genes") %in% colnames(ms)))
  expect_equal(sum(ms$n_genes), ncol(d))

  gl <- module_gene_list(wg)
  expect_false("grey" %in% names(gl))

  mod1 <- setdiff(unique(wg$modules), "grey")[1]
  hubs <- hub_genes(wg, mod1, n = 5)
  expect_true(all(c("gene", "kME") %in% colnames(hubs)))
  expect_false(is.unsorted(rev(hubs$kME)))     # kME descending
})

test_that("module_trait_cor returns aligned correlation, p and BH-adjusted p", {
  skip_if_not_installed("WGCNA")
  d <- wgcna_datexpr(sim_norm(), n_genes = 150)
  wg <- suppressWarnings(run_wgcna(d, power = 6, min_module_size = 20))
  meta <- data.frame(sample = rownames(d),
                     grp = rep(c("A", "B"), length.out = nrow(d)))
  tr <- build_traits(meta, rownames(d))
  mt <- module_trait_cor(wg$MEs, tr)
  expect_equal(dim(mt$cor), dim(mt$p))
  expect_equal(dim(mt$padj), dim(mt$p))          # FDR matrix present
  expect_true(all(mt$padj >= mt$p - 1e-9))        # BH >= raw p
  expect_equal(nrow(mt$cor), ncol(wg$MEs))
  expect_equal(mt$n, nrow(d))
  expect_true(all(mt$cor >= -1 & mt$cor <= 1))
})

test_that("wgcna_default_power follows WGCNA's sample-size recommendations", {
  expect_equal(wgcna_default_power(12, "signed"), 18)
  expect_equal(wgcna_default_power(25, "signed"), 16)
  expect_equal(wgcna_default_power(35, "signed"), 14)
  expect_equal(wgcna_default_power(50, "signed"), 12)
  expect_equal(wgcna_default_power(12, "unsigned"), 9)
})
