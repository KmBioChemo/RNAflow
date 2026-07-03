demo_path <- function(f) system.file("extdata", f, package = "RNAflow")

test_that("bundled demo datasets load and validate", {
  sets <- list(
    c("demo_airway_counts.csv", "demo_airway_metadata.csv"),
    c("demo_pickrell_counts.csv", "demo_pickrell_metadata.csv")
  )
  for (s in sets) {
    fc <- demo_path(s[1]); fm <- demo_path(s[2])
    skip_if(!file.exists(fc) || !file.exists(fm))
    counts <- read_counts(fc)
    meta   <- read_metadata(fm, counts_samples = colnames(counts))
    expect_true(is.matrix(counts))
    expect_gt(nrow(counts), 100)
    expect_identical(sort(colnames(counts)), sort(meta[[1]]))
    expect_silent(validate_counts(counts, strict = TRUE))
  }
})

test_that("airway demo has the expected real-data structure", {
  fc <- demo_path("demo_airway_counts.csv")
  fm <- demo_path("demo_airway_metadata.csv")
  skip_if(!file.exists(fc) || !file.exists(fm))
  counts <- read_counts(fc)
  meta   <- read_metadata(fm, counts_samples = colnames(counts))
  expect_equal(ncol(counts), 8)
  expect_setequal(unique(meta$condition), c("Control", "Dex"))
  expect_equal(length(unique(meta$cell)), 4)          # 4 cell lines (covariate)
  # gene IDs are symbols, so a known dexamethasone target is present
  expect_true("DUSP1" %in% rownames(counts))
})

test_that("pickrell demo has the expected real-data structure", {
  fc <- demo_path("demo_pickrell_counts.csv")
  fm <- demo_path("demo_pickrell_metadata.csv")
  skip_if(!file.exists(fc) || !file.exists(fm))
  counts <- read_counts(fc)
  meta   <- read_metadata(fm, counts_samples = colnames(counts))
  expect_equal(ncol(counts), 30)                       # balanced 15 F / 15 M
  expect_setequal(unique(meta$sex), c("female", "male"))
  # gene IDs are Ensembl (complements the symbol-based airway demo)
  expect_true(all(grepl("^ENSG", utils::head(rownames(counts), 20))))
})
