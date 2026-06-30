test_that("is_hex6 accepts valid hex codes", {
  expect_true(is_hex6("#FFFFFF"))
  expect_true(is_hex6("#000000"))
  expect_true(is_hex6("#1D9E75"))
  expect_true(is_hex6("#ff00aa"))  # lowercase
  expect_true(is_hex6(" #FF00AA "))  # whitespace
})

test_that("is_hex6 rejects invalid input", {
  expect_false(is_hex6("FFFFFF"))      # missing #
  expect_false(is_hex6("#FFF"))         # 3-digit
  expect_false(is_hex6("#GGGGGG"))      # non-hex chars
  expect_false(is_hex6("#FFFFFFF"))     # too long
  expect_false(is_hex6(NA))
  expect_false(is_hex6(NULL))
  expect_false(is_hex6(""))
  expect_false(is_hex6(c("#FFFFFF", "#000000")))  # vector of length > 1
})

test_that("safe_col returns valid hex or fallback", {
  expect_equal(safe_col("#1D9E75"), "#1D9E75")
  expect_equal(safe_col("invalid", "#FF0000"), "#FF0000")
  expect_equal(safe_col(NA, "#00FF00"), "#00FF00")
  expect_equal(safe_col(NULL, "#0000FF"), "#0000FF")
  # Both invalid: ultimate fallback
  expect_equal(safe_col("garbage", "also garbage"), "#888888")
})

test_that("make_palette returns n colors", {
  expect_length(make_palette("RdBu", 100), 100)
  expect_length(make_palette("viridis", 50), 50)
  expect_length(make_palette("Blues", 7), 7)
})

test_that("make_palette returns valid hex codes", {
  pal <- make_palette("RdBu", 10)
  expect_true(all(vapply(pal, is_hex6, logical(1))))
})

test_that("build_cols produces n valid colors even with junk input", {
  cols <- build_cols(3, c("Control", "Treatment", "Recovery"), list())
  expect_length(cols, 3)
  expect_true(all(vapply(cols, is_hex6, logical(1))))

  # With user override
  user <- list(grp_col_Control = "#FF0000")
  cols2 <- build_cols(3, c("Control", "Treatment", "Recovery"), user)
  expect_equal(cols2[1], "#FF0000")
  expect_true(all(vapply(cols2, is_hex6, logical(1))))

  # Junk user input falls back safely
  user_bad <- list(grp_col_Control = "definitely not a color")
  cols3 <- build_cols(3, c("Control", "Treatment", "Recovery"), user_bad)
  expect_length(cols3, 3)
  expect_true(all(vapply(cols3, is_hex6, logical(1))))
})

test_that("SAFE8 fallback palette is always valid", {
  expect_length(SAFE8, 8)
  expect_true(all(vapply(SAFE8, is_hex6, logical(1))))
})
