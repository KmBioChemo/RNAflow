test_that("empty_project returns canonical structure", {
  p <- empty_project("test")
  expect_equal(p$name, "test")
  expect_true(all(c("counts", "metadata", "de_results", "de_params",
                    "figures", "enrichment", "wgcna", "notes") %in% names(p)))
  expect_s3_class(p$created_at, "POSIXct")
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
