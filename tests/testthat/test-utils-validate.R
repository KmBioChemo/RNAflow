test_that("validate_counts accepts valid integer matrix", {
  m <- matrix(as.integer(rpois(60, 100)), nrow = 10, ncol = 6,
              dimnames = list(paste0("g", 1:10), paste0("s", 1:6)))
  expect_silent(validate_counts(m))
})

test_that("validate_counts rejects NULL / empty / wrong type", {
  expect_error(validate_counts(NULL), "NULL")
  expect_error(validate_counts(matrix(numeric(0), 0, 0)), "empty")
  expect_error(validate_counts("not a matrix"), "matrix or data.frame")
})

test_that("validate_counts requires rownames", {
  m <- matrix(1:12, nrow = 3, dimnames = list(NULL, paste0("s", 1:4)))
  expect_error(validate_counts(m), "rownames")
})

test_that("validate_counts requires colnames", {
  m <- matrix(1:12, nrow = 3, dimnames = list(paste0("g", 1:3), NULL))
  expect_error(validate_counts(m), "column names")
})

test_that("validate_counts catches duplicate gene IDs", {
  m <- matrix(1:12, nrow = 3, dimnames = list(c("g1", "g1", "g3"),
                                              paste0("s", 1:4)))
  expect_error(validate_counts(m), "[Dd]uplicat")
})

test_that("validate_counts catches negative values", {
  m <- matrix(c(1:11, -1), nrow = 3, dimnames = list(paste0("g", 1:3),
                                                     paste0("s", 1:4)))
  expect_error(validate_counts(m), "negative")
})

test_that("validate_counts catches non-integer in strict mode", {
  m <- matrix(c(1:11, 1.5), nrow = 3, dimnames = list(paste0("g", 1:3),
                                                      paste0("s", 1:4)))
  expect_error(validate_counts(m, strict = TRUE), "non-integer")
  # But allowed in non-strict
  expect_silent(validate_counts(m, strict = FALSE))
})

test_that("validate_counts catches NAs", {
  m <- matrix(c(1:11, NA), nrow = 3, dimnames = list(paste0("g", 1:3),
                                                     paste0("s", 1:4)))
  expect_error(validate_counts(m), "NA")
})

test_that("validate_metadata accepts a clean data.frame", {
  d <- data.frame(sample = paste0("s", 1:6),
                  condition = rep(c("Ctrl", "Trt"), each = 3),
                  stringsAsFactors = FALSE)
  expect_silent(validate_metadata(d))
})

test_that("validate_metadata rejects single-column input", {
  d <- data.frame(sample = paste0("s", 1:6))
  expect_error(validate_metadata(d), "at least 2 columns")
})

test_that("validate_metadata catches duplicate sample IDs", {
  d <- data.frame(sample = c("s1", "s1", "s2"),
                  condition = c("A", "B", "C"))
  expect_error(validate_metadata(d), "[Dd]uplicat")
})

test_that("validate_metadata flags mismatch with counts samples", {
  d <- data.frame(sample = paste0("s", 1:6),
                  condition = rep(c("Ctrl", "Trt"), each = 3))
  # All-mismatch fails outright
  expect_error(
    validate_metadata(d, counts_samples = paste0("xx", 1:6)),
    "in common"
  )
  # Partial mismatch warns
  expect_warning(
    validate_metadata(d, counts_samples = c(paste0("s", 1:4), "extra1", "extra2")),
    "not found"
  )
})

test_that("validate_de_results requires standard columns", {
  d <- data.frame(gene = "a", log2FoldChange = 1)  # missing padj
  expect_error(validate_de_results(d), "padj")
})

test_that("validate_de_results coerces numeric columns", {
  d <- data.frame(gene = c("a", "b"),
                  log2FoldChange = c("1.5", "-2"),
                  padj = c("0.01", "0.001"),
                  stringsAsFactors = FALSE)
  expect_silent(validate_de_results(d))
})
