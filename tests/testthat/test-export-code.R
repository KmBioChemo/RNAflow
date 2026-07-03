make_project <- function(n = 2) {
  p <- empty_project("demo")
  p$organism <- "mouse"
  res <- data.frame(gene = c("a", "b"), log2FoldChange = c(1, -1),
                    padj = c(0.01, 0.02))
  if (n >= 1)
    p$contrasts <- contrast_store_upsert(
      p$contrasts, "group: WT_LPS vs WT_Veh", res,
      list(design_var = "group", treated = "WT_LPS", reference = "WT_Veh",
           shrink = FALSE, min_count = 10, alpha = 0.05))
  if (n >= 2)
    p$contrasts <- contrast_store_upsert(
      p$contrasts, "group: KO_LPS vs KO_Veh", res,
      list(design_var = "group", treated = "KO_LPS", reference = "KO_Veh",
           shrink = TRUE, min_count = 10, alpha = 0.05))
  p
}

test_that("generated script is syntactically valid R", {
  code <- generate_r_script(make_project(2))
  expect_type(code, "character")
  expect_length(code, 1)
  expect_silent(parse(text = code))          # must parse
})

test_that("script reflects the stored contrasts and organism", {
  code <- generate_r_script(make_project(2), counts_path = "my_counts.csv")
  expect_match(code, "library(RNAflow)", fixed = TRUE)
  expect_match(code, 'read_counts("my_counts.csv")', fixed = TRUE)
  expect_match(code, 'organism <- "mouse"', fixed = TRUE)
  expect_match(code, 'contrast = c("group", "WT_LPS", "WT_Veh")', fixed = TRUE)
  expect_match(code, 'contrast = c("group", "KO_LPS", "KO_Veh")', fixed = TRUE)
  # multi-contrast block present with 2 contrasts
  expect_match(code, "fig_upset(contrast_sig_sets", fixed = TRUE)
  # downstream analyses templated
  expect_match(code, "run_gsea", fixed = TRUE)
  expect_match(code, "run_wgcna", fixed = TRUE)
  expect_match(code, "sessionInfo()", fixed = TRUE)
})

test_that("empty project still yields a valid template script", {
  code <- generate_r_script(empty_project("empty"))
  expect_silent(parse(text = code))
  expect_match(code, "No contrasts were saved", fixed = TRUE)
  expect_no_match(code, "fig_upset", fixed = TRUE)   # no compare block
})

test_that("generated timestamp is stamped when provided", {
  code <- generate_r_script(make_project(1), generated = "2026-07-01 10:00")
  expect_match(code, "on 2026-07-01 10:00", fixed = TRUE)
})
