# Linked-explorer data prep (pure) + figure smoke test.

de_lk <- data.frame(
  gene = paste0("g", 1:6),
  log2FoldChange = c(3, -2.5, 0.2, 1.5, -1.4, NA),
  stat = c(8, -7, 0.5, 3, -3, 0),
  pvalue = c(1e-9, 1e-8, 0.4, 1e-3, 1e-3, 0.5),
  padj = c(1e-8, 1e-7, 0.5, 2e-3, 3e-3, NA),
  stringsAsFactors = FALSE)

test_that("linked_volcano_df categorizes significance and drops NAs", {
  df <- linked_volcano_df(de_lk, padj_thr = 0.05, lfc_thr = 1)
  expect_true(all(c("gene", "log2FoldChange", "negLog10P", "padj",
                    "significance") %in% colnames(df)))
  expect_equal(nrow(df), 5)                       # g6 (NA padj) dropped
  sig <- stats::setNames(df$significance, df$gene)
  expect_equal(sig[["g1"]], "Up")
  expect_equal(sig[["g2"]], "Down")
  expect_equal(sig[["g3"]], "NS")                 # |lfc| < 1
  expect_true(all(df$negLog10P >= 0))
})

test_that("linked_volcano_df respects custom thresholds", {
  df <- linked_volcano_df(de_lk, padj_thr = 0.05, lfc_thr = 2)
  sig <- stats::setNames(df$significance, df$gene)
  expect_equal(sig[["g4"]], "NS")                 # lfc 1.5 < 2 now
  expect_equal(sig[["g1"]], "Up")                 # lfc 3 still up
})

test_that("fig_linked_volcano returns a plotly widget from SharedData", {
  skip_if_not_installed("crosstalk")
  skip_if_not_installed("plotly")
  df <- linked_volcano_df(de_lk)
  sd <- crosstalk::SharedData$new(df, key = ~gene)
  p <- fig_linked_volcano(sd)
  expect_s3_class(p, "plotly")
})
