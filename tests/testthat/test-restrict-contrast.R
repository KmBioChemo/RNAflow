meta6 <- data.frame(
  sample = paste0("s", 1:12),
  group  = rep(c("A", "B", "C"), each = 4),
  batch  = rep(c("X", "Y"), 6),
  stringsAsFactors = FALSE)
counts6 <- matrix(1L, nrow = 20, ncol = 12,
                  dimnames = list(paste0("g", 1:20), paste0("s", 1:12)))

test_that("restrict_to_contrast keeps only the two contrast groups", {
  params <- list(design_var = "group", treated = "A", reference = "B")
  out <- restrict_to_contrast(counts6, meta6, params)
  expect_equal(ncol(out$counts), 8)                         # A (4) + B (4)
  expect_setequal(colnames(out$counts), paste0("s", 1:8))
  expect_equal(nrow(out$metadata), 8)
  expect_setequal(unique(out$metadata$group), c("A", "B"))
})

test_that("restrict_to_contrast returns inputs unchanged when not applicable", {
  # missing params
  expect_equal(ncol(restrict_to_contrast(counts6, meta6, NULL)$counts), 12)
  # design variable absent from metadata
  bad <- list(design_var = "nope", treated = "A", reference = "B")
  expect_equal(ncol(restrict_to_contrast(counts6, meta6, bad)$counts), 12)
  # uploaded-style params (no design_var)
  up <- list(source = "upload")
  expect_equal(ncol(restrict_to_contrast(counts6, meta6, up)$counts), 12)
})
