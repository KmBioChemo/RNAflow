# Pure prompt-building / helpers for the AI interpretation layer.
# No network access -- call_claude() is only exercised for its guard clauses.

de_ai <- data.frame(
  gene = c("AAA", "BBB", "CCC", "DDD", "EEE", "FFF"),
  log2FoldChange = c(3.0, 2.0, 0.2, -2.5, -1.8, -0.1),
  stat  = c(8, 6, 0.5, -7, -5, -0.3),
  pvalue = c(1e-9, 1e-6, 0.4, 1e-8, 1e-5, 0.6),
  padj  = c(1e-8, 1e-5, 0.5, 1e-7, 1e-4, 0.7),
  stringsAsFactors = FALSE)

test_that("summarize_de_for_ai counts up/down and lists the right genes", {
  s <- summarize_de_for_ai(de_ai, top_n = 30, padj_thr = 0.05, lfc_thr = 1)
  expect_true(grepl("2 up, 2 down", s))          # AAA,BBB up ; DDD,EEE down
  expect_true(grepl("AAA", s) && grepl("DDD", s))
  expect_false(grepl("CCC", s))                  # not significant
})

test_that("summarize_de_for_ai respects top_n per direction", {
  s <- summarize_de_for_ai(de_ai, top_n = 1)
  # Only the single strongest up gene (AAA) is listed, not BBB
  expect_true(grepl("AAA", s))
  expect_false(grepl("BBB", s))
})

test_that("summarize_enrich_for_ai handles NULL, empty, GSEA and ORA", {
  expect_null(summarize_enrich_for_ai(NULL))
  expect_null(summarize_enrich_for_ai(list(method = "gsea",
                                           table = data.frame())))

  gsea <- list(method = "gsea", table = data.frame(
    pathway = c("HALLMARK_INFLAMMATION", "HALLMARK_HYPOXIA"),
    NES = c(2.3, -1.9), padj = c(1e-5, 1e-3), stringsAsFactors = FALSE))
  g <- summarize_enrich_for_ai(gsea)
  expect_true(grepl("HALLMARK_INFLAMMATION", g))
  expect_true(grepl("NES", g))

  ora <- list(method = "ora", table = data.frame(
    Description = c("immune response", "cell cycle"),
    Count = c(12, 8), padj = c(1e-4, 1e-2), stringsAsFactors = FALSE))
  o <- summarize_enrich_for_ai(ora)
  expect_true(grepl("immune response", o))
  expect_true(grepl("genes=12", o))
})

test_that("build_interpret_prompt embeds contrast, genes and enrichment", {
  gsea <- list(method = "gsea", table = data.frame(
    pathway = "HALLMARK_INFLAMMATION", NES = 2.3, padj = 1e-5,
    stringsAsFactors = FALSE))
  p <- build_interpret_prompt(
    de_ai, enrich = gsea, organism = "human",
    contrast_params = list(treated = "dex", reference = "ctrl",
                           design_var = "condition"))
  expect_true(grepl("Homo sapiens", p))          # organism taxon
  expect_true(grepl("dex vs ctrl", p))           # contrast direction
  expect_true(grepl("AAA", p))                   # a top gene
  expect_true(grepl("HALLMARK_INFLAMMATION", p)) # enrichment term
  expect_true(grepl("Summary", p))               # requested section
})

test_that("build_interpret_prompt works without enrichment", {
  p <- build_interpret_prompt(de_ai, enrich = NULL, organism = "mouse")
  expect_true(grepl("No functional enrichment", p))
  expect_true(grepl("Mus musculus", p))
})

test_that("estimate_cost computes and degrades gracefully", {
  # 1M in + 1M out on Opus 4.8 = $5 + $25 = $30
  expect_equal(estimate_cost(1e6, 1e6, "claude-opus-4-8"), 30)
  expect_true(is.na(estimate_cost(NA, 100, "claude-opus-4-8")))
  expect_true(is.na(estimate_cost(100, 100, "unknown-model")))
})

test_that("call_claude refuses without an API key", {
  expect_error(call_claude("hi", api_key = ""), "API key")
})
