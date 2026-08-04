test_that("read_counts collapses duplicate gene IDs by summing (default)", {
  tmp <- tempfile(fileext = ".csv")
  df <- data.frame(gene = c("A", "B", "A", "C"),
                   s1 = c(1L, 2L, 4L, 3L),
                   s2 = c(10L, 20L, 40L, 30L))
  utils::write.csv(df, tmp, row.names = FALSE)

  m <- read_counts(tmp)
  expect_equal(nrow(m), 3L)
  expect_false(anyDuplicated(rownames(m)) > 0)
  expect_equal(unname(m["A", "s1"]), 5)    # 1 + 4
  expect_equal(unname(m["A", "s2"]), 50)   # 10 + 40
  expect_equal(attr(m, "n_collapsed"), 1L)
  # summed counts stay integer-valued
  expect_true(all(m == round(m)))
})

test_that("read_counts can keep the most-expressed duplicate", {
  tmp <- tempfile(fileext = ".csv")
  df <- data.frame(gene = c("A", "A", "B"),
                   s1 = c(1L, 100L, 5L),
                   s2 = c(1L, 100L, 5L))
  utils::write.csv(df, tmp, row.names = FALSE)

  m <- read_counts(tmp, duplicate_action = "max")
  expect_equal(nrow(m), 2L)
  expect_equal(unname(m["A", "s1"]), 100)  # kept the higher-total row
})

test_that("read_counts still rejects duplicates when asked", {
  tmp <- tempfile(fileext = ".csv")
  df <- data.frame(gene = c("A", "A"), s1 = c(1L, 2L), s2 = c(3L, 4L))
  utils::write.csv(df, tmp, row.names = FALSE)

  expect_error(read_counts(tmp, duplicate_action = "reject"),
               "Duplicated gene IDs")
})

test_that("unique gene IDs are unaffected and report zero collapses", {
  tmp <- tempfile(fileext = ".csv")
  df <- data.frame(gene = c("A", "B", "C"),
                   s1 = c(1L, 2L, 3L), s2 = c(4L, 5L, 6L))
  utils::write.csv(df, tmp, row.names = FALSE)

  m <- read_counts(tmp)
  expect_equal(nrow(m), 3L)
  expect_equal(attr(m, "n_collapsed"), 0L)
})
