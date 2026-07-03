make_methods_project <- function() {
  assemble_project(
    "t", "human", NULL, NULL,
    contrasts = contrast_store_upsert(
      list(), "condition: Dex vs Control",
      data.frame(gene = "g", log2FoldChange = 1, padj = 0.01),
      list(design_var = "condition", treated = "Dex", reference = "Control",
           covariates = "cell", shrink = TRUE, shrink_used = "apeglm",
           min_count = 10, alpha = 0.05)),
    settings = list(
      enrichment = list(method = "gsea", collection = "H", rank_by = "stat"),
      wgcna = list(n_genes = 3000, network_type = "signed", power = 12,
                   min_module_size = 30, merge_cut_height = 0.25)))
}

test_that("generate_methods_text summarizes DE, enrichment and WGCNA", {
  txt <- generate_methods_text(make_methods_project())
  expect_type(txt, "character"); expect_length(txt, 1)
  expect_match(txt, "DESeq2")
  expect_match(txt, "adjusting for cell", fixed = TRUE)
  expect_match(txt, "Dex vs Control", fixed = TRUE)
  expect_match(txt, "apeglm", fixed = TRUE)
  expect_match(txt, "unshrunken model", fixed = TRUE)  # inference separation
  expect_match(txt, "fgsea", fixed = TRUE)
  expect_match(txt, "Wald statistic", fixed = TRUE)
  expect_match(txt, "WGCNA", fixed = TRUE)
  expect_match(txt, "power of 12", fixed = TRUE)
  expect_match(txt, "RNAflow", fixed = TRUE)
})

test_that("generate_methods_text handles an ORA run and an empty project", {
  p <- make_methods_project()
  p$enrichment <- list(method = "ora", db = "GO", ont = "BP", padj = 0.05, lfc = 1)
  txt <- generate_methods_text(p)
  expect_match(txt, "Over-representation analysis", fixed = TRUE)
  expect_match(txt, "clusterProfiler", fixed = TRUE)
  expect_match(txt, "BP ontology", fixed = TRUE)

  empty <- generate_methods_text(empty_project("x"))
  expect_match(empty, "No analysis has been run", fixed = TRUE)
})
