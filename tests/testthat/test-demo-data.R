demo_path <- function(f) system.file("extdata", f, package = "RNAflow")

test_that("bundled demo datasets load and validate", {
  sets <- list(
    c("demo_airway_counts.csv", "demo_airway_metadata.csv"),
    c("demo_tcga_counts.csv", "demo_tcga_metadata.csv")
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

test_that("TCGA demo has the expected complex multi-group structure", {
  fc <- demo_path("demo_tcga_counts.csv")
  fm <- demo_path("demo_tcga_metadata.csv")
  skip_if(!file.exists(fc) || !file.exists(fm))
  counts <- read_counts(fc)
  meta   <- read_metadata(fm, counts_samples = colnames(counts))
  expect_equal(ncol(counts), 120)                      # 8 types x 15 tumors
  expect_equal(length(unique(meta$cancer_type)), 8)    # eight cancer types
  expect_setequal(unique(meta$cancer_type),
                  c("BRCA", "LUAD", "KIRC", "LGG", "THCA", "PRAD", "COAD", "SKCM"))
  # gene IDs are symbols, so a canonical cancer gene is present
  expect_true("TP53" %in% rownames(counts))
})
