make_report_project <- function() {
  set.seed(5)
  mkres <- function(n = 300) data.frame(
    gene = paste0("g", seq_len(n)), baseMean = rgamma(n, 2, 0.1),
    log2FoldChange = rnorm(n, 0, 2), lfcSE = 0.2, stat = rnorm(n),
    pvalue = runif(n), padj = pmin(runif(n) * 0.4, 1),
    stringsAsFactors = FALSE)
  p <- empty_project("report demo")
  p$organism <- "mouse"
  p$counts   <- matrix(1L, nrow = 300, ncol = 6,
                       dimnames = list(paste0("g", 1:300), paste0("s", 1:6)))
  p$metadata <- data.frame(sample = paste0("s", 1:6),
                           group = rep(c("A", "B"), each = 3))
  p$contrasts <- contrast_store_upsert(p$contrasts, "group: B vs A", mkres(),
    list(design_var = "group", treated = "B", reference = "A",
         shrink = FALSE, min_count = 10, alpha = 0.05))
  p$contrasts <- contrast_store_upsert(p$contrasts, "group: A vs B", mkres(),
    list(design_var = "group", treated = "A", reference = "B",
         shrink = FALSE, min_count = 10, alpha = 0.05))
  p
}

test_that("session_manifest lists RNAflow and installed packages", {
  m <- session_manifest()
  expect_s3_class(m, "data.frame")
  expect_true("RNAflow" %in% m$package)
  expect_true(all(nzchar(m$version)))
})

test_that("build_report_html writes a self-contained HTML file", {
  skip_if_not_installed("xfun")
  f <- tempfile(fileext = ".html")
  out <- build_report_html(make_report_project(), f, generated = "2026-07-01")
  expect_equal(out, f)
  expect_true(file.exists(f))
  expect_gt(file.size(f), 5000)              # non-trivial (embedded figures)

  html <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_match(html, "RNAflow analysis report", fixed = TRUE)
  expect_match(html, "Reproducible R script", fixed = TRUE)
  expect_match(html, "group: B vs A", fixed = TRUE)
  expect_match(html, "data:image/png;base64,", fixed = TRUE)   # embedded figure
  expect_match(html, "Cross-contrast signature", fixed = TRUE) # 2 contrasts
})

test_that("build_report_html works with an empty project", {
  f <- tempfile(fileext = ".html")
  expect_silent(build_report_html(empty_project("empty"), f))
  expect_true(file.exists(f))
})

test_that("build_report_html renders an AI interpretation when present", {
  p <- make_report_project()
  p$ai_interpretation <- list(
    text = "## Summary\nStrong **GR activation** signature.",
    model = "claude-opus-4-8")
  f <- tempfile(fileext = ".html")
  build_report_html(p, f)
  html <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_match(html, "AI interpretation", fixed = TRUE)
  expect_match(html, "GR activation", fixed = TRUE)
  expect_match(html, "claude-opus-4-8", fixed = TRUE)
})

test_that("assemble_project carries the AI interpretation from settings", {
  s <- list(ai_interpretation = list(text = "hello", model = "claude-opus-4-8"))
  p <- assemble_project("demo", organism = "human", settings = s)
  expect_equal(p$ai_interpretation$text, "hello")
})

test_that("build_report_html renders a Signatures section when recorded", {
  skip_if_not_installed("xfun")
  skip_if_not_installed("pheatmap")
  p <- make_report_project()
  set.seed(2)
  es <- matrix(rnorm(5 * 6), nrow = 5,
               dimnames = list(paste0("SET_", 1:5), paste0("s", 1:6)))
  p$signatures <- list(collection = "MSigDB Hallmark", method = "gsva",
                       organism = "mouse", group_by = "group",
                       n_top = 5, n_sets = 5, n_samples = 6,
                       generated = "2026-07-03 10:00:00", scores = es)
  f <- tempfile(fileext = ".html")
  build_report_html(p, f)
  html <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_match(html, "Signatures (GSVA", fixed = TRUE)
})
