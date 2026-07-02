# Activity inference (decoupleR). Pure helpers are tested directly; the
# network-fetch + scoring functions are guarded and only run if decoupleR is
# installed.

de_act <- data.frame(
  gene = paste0("G", 1:8),
  log2FoldChange = c(3, 2, 1, 0.5, -0.5, -1, -2, -3),
  stat = c(9, 6, 3, 1, -1, -3, -6, -9),
  pvalue = 10^(-c(9, 6, 3, 1, 1, 3, 6, 9)),
  padj = 10^(-c(8, 5, 2, 1, 1, 2, 5, 8)),
  stringsAsFactors = FALSE)

test_that("decoupler_organism maps supported organisms and rejects others", {
  expect_equal(decoupler_organism("Human"), "human")
  expect_equal(decoupler_organism("mouse"), "mouse")
  expect_error(decoupler_organism("frog"), "human")
})

test_that("activity_input builds a one-column ranked matrix", {
  m <- activity_input(de_act, by = "stat")
  expect_true(is.matrix(m))
  expect_equal(ncol(m), 1L)
  expect_equal(nrow(m), 8L)
  expect_equal(rownames(m)[1], "G1")          # highest stat first
  expect_true(m[1, 1] > m[nrow(m), 1])        # sorted decreasing
})

test_that("activity_input errors on too few genes", {
  expect_error(activity_input(de_act[1, , drop = FALSE]), "Too few")
})

test_that("run_activity validates its network argument", {
  skip_if_not_installed("decoupleR")
  expect_error(run_activity(de_act, data.frame()), "non-empty")
  expect_error(
    run_activity(de_act, data.frame(source = "T", target = "G1")),
    "weight column")
})

test_that("run_activity scores a synthetic regulon network (ulm)", {
  skip_if_not_installed("decoupleR")
  net <- data.frame(
    source = rep(c("TF_UP", "TF_DN"), each = 4),
    target = c("G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8"),
    mor = c(1, 1, 1, 1, 1, 1, 1, 1),
    stringsAsFactors = FALSE)
  out <- run_activity(de_act, net, method = "ulm", mor_col = "mor",
                      min_size = 3)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("source", "score", "p_value", "padj") %in% colnames(out)))
  # TF_UP targets the top (positive) genes -> positive activity;
  # TF_DN targets the bottom (negative) genes -> negative activity.
  sc <- stats::setNames(out$score, out$source)
  expect_gt(sc[["TF_UP"]], sc[["TF_DN"]])
})

test_that("run_activity scores a synthetic pathway network (mlm)", {
  skip_if_not_installed("decoupleR")
  # PROGENy-style network with a 'weight' column, scored with the multivariate
  # model -- no OmniPath / internet needed. The multivariate fit needs
  # non-collinear footprints, so weights vary and 4 background genes (H9-H12)
  # sit outside both pathways.
  de12 <- data.frame(
    gene = paste0("H", 1:12),
    log2FoldChange = c(3, 2.5, 2, 1.5, -1.5, -2, -2.5, -3, 0.2, -0.1, 0.1, 0),
    stat = c(9, 7, 6, 4, -4, -6, -7, -9, 0.5, -0.3, 0.4, 0),
    pvalue = 10^(-c(9, 7, 6, 4, 4, 6, 7, 9, 1, 1, 1, 1)),
    padj = 10^(-c(8, 6, 5, 3, 3, 5, 6, 8, 1, 1, 1, 1)),
    stringsAsFactors = FALSE)
  net <- data.frame(
    source = c(rep("PW_UP", 4), rep("PW_DN", 4)),
    target = c("H1", "H2", "H3", "H4", "H5", "H6", "H7", "H8"),
    weight = c(2, 1.5, 1, 0.5, 2, 1.5, 1, 0.5),
    stringsAsFactors = FALSE)
  out <- run_activity(de12, net, method = "mlm", mor_col = "weight",
                      min_size = 3)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("source", "score", "p_value", "padj") %in% colnames(out)))
  sc <- stats::setNames(out$score, out$source)
  expect_gt(sc[["PW_UP"]], sc[["PW_DN"]])   # up-footprint activated, down repressed
})

test_that("fig_activity_bar returns a ggplot for a synthetic activity table", {
  act <- data.frame(
    source = paste0("TF", 1:6),
    score = c(4, -3, 2.5, -2, 1.5, -1),
    padj = 10^(-c(6, 5, 4, 3, 2, 1)), stringsAsFactors = FALSE)
  p <- fig_activity_bar(act, n = 5)
  expect_s3_class(p, "ggplot")
  expect_error(fig_activity_bar(act[0, ]), "No activity")
})
