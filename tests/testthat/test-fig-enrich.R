fake_gsea <- function() {
  data.frame(
    pathway = c("HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_E2F_TARGETS",
                "HALLMARK_OXIDATIVE_PHOSPHORYLATION"),
    pval = c(1e-20, 1e-5, 1e-3),
    padj = c(1e-18, 1e-4, 1e-2),
    NES = c(2.9, -1.8, 1.2),
    ES = c(0.7, -0.5, 0.4),
    size = c(200L, 180L, 150L),
    stringsAsFactors = FALSE
  )
}

fake_ora <- function() {
  data.frame(
    ID = c("GO:1", "GO:2"),
    Description = c("response to lipopolysaccharide", "response to virus"),
    GeneRatio = c("66/402", "61/402"),
    BgRatio = c("300/20000", "280/20000"),
    pvalue = c(1e-30, 1e-25),
    padj = c(1e-28, 1e-23),
    qvalue = c(1e-28, 1e-23),
    Count = c(66L, 61L),
    geneID = c("a/b/c", "d/e/f"),
    stringsAsFactors = FALSE
  )
}

test_that("enrich_plot_df normalizes GSEA and ORA tables", {
  g <- enrich_plot_df(fake_gsea())
  expect_true(all(c("term", "score", "padj", "direction", "is_gsea") %in% names(g)))
  expect_true(all(g$is_gsea))
  expect_setequal(unique(g$direction), c("Up", "Down"))

  o <- enrich_plot_df(fake_ora())
  expect_false(o$is_gsea[1])
  expect_true(all(o$score > 0 & o$score < 1))   # parsed GeneRatio
})

test_that("enrich_plot_df errors on empty input", {
  expect_error(enrich_plot_df(fake_gsea()[0, ]), "No enrichment results")
})

test_that("clean_term strips HALLMARK_ and underscores", {
  expect_equal(clean_term("HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
               "TNFA SIGNALING VIA NFKB")
})

test_that("fig_enrich_dot returns a ggplot for GSEA and ORA", {
  expect_s3_class(fig_enrich_dot(fake_gsea()), "ggplot")
  expect_s3_class(fig_enrich_dot(fake_ora()), "ggplot")
  expect_s3_class(fig_enrich_dot(fake_gsea(), mode = "publication"), "ggplot")
})

test_that("fig_enrich_bar returns a ggplot for GSEA and ORA", {
  expect_s3_class(fig_enrich_bar(fake_gsea()), "ggplot")
  expect_s3_class(fig_enrich_bar(fake_ora()), "ggplot")
})

test_that("fig_gsea_curve builds from DE results and a gene set", {
  skip_if_not_installed("fgsea")
  set.seed(3)
  n <- 800
  res <- data.frame(gene = paste0("G", seq_len(n)),
                    log2FoldChange = rnorm(n), padj = runif(n),
                    pvalue = runif(n), stat = rnorm(n))
  res <- res[order(res$stat, decreasing = TRUE), ]
  pg <- head(res$gene, 50)
  expect_s3_class(fig_gsea_curve(res, pg, title = "HALLMARK_TEST"), "ggplot")
})
