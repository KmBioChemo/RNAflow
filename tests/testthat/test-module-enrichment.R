fake_combined <- function() {
  do.call(rbind, lapply(c("blue", "brown", "yellow"), function(m) {
    data.frame(
      ID = paste0("GO:", m, 1:4),
      Description = paste(m, "term", 1:4),
      GeneRatio = "20/300", BgRatio = "100/20000",
      pvalue = 10^(-(10:7)), padj = 10^(-(9:6)), qvalue = 10^(-(9:6)),
      Count = c(40L, 30L, 25L, 20L), geneID = "a/b/c",
      module = m, stringsAsFactors = FALSE)
  }))
}

test_that("fig_module_enrichment returns a ggplot", {
  expect_s3_class(fig_module_enrichment(fake_combined()), "ggplot")
  expect_s3_class(fig_module_enrichment(fake_combined(), mode = "publication"),
                  "ggplot")
})

test_that("fig_module_enrichment respects max_terms", {
  p <- fig_module_enrichment(fake_combined(), max_terms = 5)
  # data behind the plot keeps at most 5 distinct terms
  expect_lte(length(unique(p$data$Description)), 5)
})

test_that("fig_module_enrichment validates its input", {
  expect_error(fig_module_enrichment(data.frame(x = 1)), "enrich_modules")
  expect_error(fig_module_enrichment(fake_combined()[0, ]), "enrich_modules")
})

test_that("enrich_modules maps each module back to real biology", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("msigdbr")
  skip_on_cran()
  # Two real Hallmark gene sets assigned to two modules. enrich_modules should
  # map each module (symbols -> ENTREZ) and recover coherent GO terms. A
  # hand-built WGCNA object keeps this deterministic -- WGCNA's own module
  # detection on a small simulated matrix is realization-dependent and not what
  # this test is about (that is covered in test-analysis-wgcna.R).
  sets <- suppressMessages(get_gene_sets("human", collection = "H"))
  g1 <- utils::head(unique(sets[["HALLMARK_INTERFERON_GAMMA_RESPONSE"]]), 45)
  g2 <- utils::head(unique(sets[["HALLMARK_E2F_TARGETS"]]), 45)
  skip_if(length(g1) < 30 || length(g2) < 30, "Hallmark sets unavailable")
  genes <- c(g1, g2)

  wg <- list(
    modules = stats::setNames(
      c(rep("blue", length(g1)), rep("brown", length(g2))), genes),
    datExpr = matrix(0, nrow = 2, ncol = length(genes),
                     dimnames = list(c("s1", "s2"), genes)))

  comb <- suppressMessages(suppressWarnings(
    enrich_modules(wg, "human", db = "GO", ont = "BP", n_per = 3)))
  expect_true("module" %in% colnames(comb))
  expect_gt(nrow(comb), 0)
  expect_setequal(unique(comb$module), c("blue", "brown"))   # both modules enrich
  expect_s3_class(fig_module_enrichment(comb), "ggplot")
})
