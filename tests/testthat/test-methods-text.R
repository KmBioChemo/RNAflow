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

test_that("generate_methods_text flags contrasts with differing settings", {
  s <- contrast_store_upsert(
    list(), "c1", data.frame(gene = "g", log2FoldChange = 1, padj = 0.01),
    list(design_var = "condition", treated = "A", reference = "C",
         shrink_used = "apeglm", min_count = 10))
  s <- contrast_store_upsert(
    s, "c2", data.frame(gene = "g", log2FoldChange = 1, padj = 0.01),
    list(design_var = "condition", treated = "B", reference = "C",
         shrink_used = "none", min_count = 50))
  txt <- generate_methods_text(assemble_project("t", "human", contrasts = s))
  expect_match(txt, "varied between contrasts", fixed = TRUE)
  expect_false(grepl("count below", txt))   # no single shared threshold claimed
})

test_that("generate_methods_text records normalization and a version caveat", {
  p <- assemble_project(
    "t", "human",
    settings = list(
      wgcna = list(n_genes = 3000, network_type = "signed", power = 6),
      normalization_method = "log2(counts+1) [VST fallback]"))
  txt <- generate_methods_text(p)
  expect_match(txt, "log2(counts+1) [VST fallback]-normalized", fixed = TRUE)
  expect_match(txt, "reflect the environment at the time", fixed = TRUE)
})
