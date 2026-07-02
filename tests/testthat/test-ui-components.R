test_that("ui_banner renders the right variant class and content", {
  b <- ui_banner("hello world", type = "warning")
  expect_s3_class(b, "shiny.tag")
  html <- as.character(b)
  expect_match(html, "rnaflow-banner", fixed = TRUE)
  expect_match(html, "rf-warning", fixed = TRUE)
  expect_match(html, "hello world", fixed = TRUE)

  # default type is info (no variant modifier), icon can be suppressed
  info <- as.character(ui_banner("x", icon = FALSE))
  expect_match(info, "rnaflow-banner", fixed = TRUE)
  expect_false(grepl("rf-warning|rf-danger|rf-success", info))
})

test_that("ui_empty_state includes title, message and icon", {
  e <- ui_empty_state("No data yet", "Upload a counts matrix to begin.",
                      icon = "inbox")
  html <- as.character(e)
  expect_match(html, "rnaflow-empty", fixed = TRUE)
  expect_match(html, "No data yet", fixed = TRUE)
  expect_match(html, "Upload a counts matrix", fixed = TRUE)
})

test_that("ui_page_header and ui_stat_tile produce expected markup", {
  h <- as.character(ui_page_header("Volcano", "Differential expression"))
  expect_match(h, "rnaflow-page-header", fixed = TRUE)
  expect_match(h, "Differential expression", fixed = TRUE)

  t <- as.character(ui_stat_tile(42, "up", bg = "#1D9E75"))
  expect_match(t, "rnaflow-tile", fixed = TRUE)
  expect_match(t, "42", fixed = TRUE)
  expect_match(t, "#1D9E75", fixed = TRUE)
})
