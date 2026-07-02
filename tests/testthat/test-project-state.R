test_that("empty_project returns canonical structure", {
  p <- empty_project("test")
  expect_equal(p$name, "test")
  expect_true(all(c("counts", "metadata", "de_results", "de_params",
                    "figures", "enrichment", "wgcna", "activity",
                    "ai_interpretation", "notes") %in% names(p)))
  expect_s3_class(p$created_at, "POSIXct")
})

test_that("assemble_project carries activity and AI settings", {
  s <- list(
    activity = list(type = "pathway", table = data.frame(source = "MAPK")),
    ai_interpretation = list(text = "hi", model = "claude-opus-4-8"))
  p <- assemble_project("demo", organism = "human", settings = s)
  expect_equal(p$activity$type, "pathway")
  expect_equal(p$ai_interpretation$model, "claude-opus-4-8")
  # No settings -> empty defaults, not an error
  p2 <- assemble_project("demo2", organism = "human")
  expect_equal(p2$activity, list())
  expect_null(p2$ai_interpretation)
})

test_that("old projects without new slots still load and render", {
  # A project saved before the activity / ai_interpretation slots existed.
  old <- empty_project("legacy")
  old$activity <- NULL; old$ai_interpretation <- NULL
  old$organism <- "mouse"
  tf <- tempfile()
  q <- load_project(save_project(old, tf))
  expect_equal(q$name, "legacy")
  # build_report_html must not error on the missing slots
  f <- tempfile(fileext = ".html")
  expect_silent(build_report_html(q, f))
  expect_true(file.exists(f))
})

test_that("save_project and load_project roundtrip", {
  skip_if_not_installed("RNAflow")
  p <- empty_project("roundtrip")
  p$notes <- c("first note", "second note")
  p$organism <- "mouse"

  tf <- tempfile()
  path <- save_project(p, tf)
  expect_true(file.exists(path))
  expect_true(grepl("\\.rnaflow\\.rds$", path))

  q <- load_project(path)
  expect_equal(q$name, "roundtrip")
  expect_equal(q$organism, "mouse")
  expect_equal(q$notes, c("first note", "second note"))
})

test_that("load_project rejects non-RNAflow files", {
  tf <- tempfile()
  saveRDS(list(some = "random object"), tf)
  expect_error(load_project(tf), "valid RNAflow")
})

test_that("contrast_store_upsert adds, updates in place, and keeps order", {
  res1 <- data.frame(gene = "g1", log2FoldChange = 1, padj = 0.01)
  res2 <- data.frame(gene = "g2", log2FoldChange = 2, padj = 0.02)

  s <- contrast_store_upsert(list(), "A", res1)
  s <- contrast_store_upsert(s, "B", res2)
  expect_named(s, c("A", "B"))
  expect_equal(s$A$results$gene, "g1")

  # Re-adding "A" updates in place without changing position
  res1b <- data.frame(gene = "g1b", log2FoldChange = 9, padj = 0.001)
  s <- contrast_store_upsert(s, "A", res1b)
  expect_named(s, c("A", "B"))
  expect_equal(s$A$results$gene, "g1b")

  expect_error(contrast_store_upsert(list(), "", res1), "non-empty")
})

test_that("contrast_store_results strips to a named list of data.frames", {
  res1 <- data.frame(gene = "g1", log2FoldChange = 1, padj = 0.01)
  s <- contrast_store_upsert(list(), "A", res1)
  out <- contrast_store_results(s)
  expect_named(out, "A")
  expect_s3_class(out$A, "data.frame")
  expect_equal(contrast_store_results(list()), list())
})

test_that("recent project caching lists newest first by project name", {
  withr::local_options(rnaflow.recent_dir = withr::local_tempdir())
  expect_equal(nrow(list_recent_projects()), 0)

  p1 <- empty_project("alpha"); p1$modified_at <- as.POSIXct("2026-01-01 10:00")
  f1 <- tempfile(fileext = ".rnaflow.rds"); saveRDS(p1, f1)
  cache_recent_project(f1, p1$name)

  p2 <- empty_project("beta"); p2$modified_at <- as.POSIXct("2026-06-01 10:00")
  f2 <- tempfile(fileext = ".rnaflow.rds"); saveRDS(p2, f2)
  cache_recent_project(f2, p2$name)

  rp <- list_recent_projects()
  expect_equal(nrow(rp), 2)
  expect_equal(rp$name[1], "beta")   # newest modified_at first
  expect_true(all(file.exists(rp$file)))
})
