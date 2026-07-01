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
