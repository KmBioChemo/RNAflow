ora_df <- function() data.frame(
  Description = paste0("term ", 1:5),
  padj  = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1),
  Count = c(20, 15, 10, 8, 5),
  geneID = c("A/B/C/D", "B/C/D/E", "C/D/E/F", "X/Y/Z", "M/N"),
  stringsAsFactors = FALSE)

gsea_df <- function() data.frame(
  pathway = paste0("p", 1:5),
  NES  = c(2.1, -1.8, 1.2, -0.9, 1.5),
  padj = c(1e-4, 1e-3, 1e-2, 2e-2, 3e-2),
  size = c(30, 25, 20, 15, 10),
  leadingEdge = I(list(c("A", "B", "C"), c("B", "C", "D"),
                       c("C", "D"), c("X"), c("M", "N"))),
  stringsAsFactors = FALSE)

test_that("fig_enrich_visnet builds a visNetwork widget for ORA and GSEA", {
  skip_if_not_installed("visNetwork")
  w1 <- fig_enrich_visnet(ora_df(), n = 5, min_similarity = 0.1)
  expect_s3_class(w1, "visNetwork")
  expect_s3_class(w1, "htmlwidget")

  w2 <- fig_enrich_visnet(gsea_df(), n = 5, min_similarity = 0.1)
  expect_s3_class(w2, "visNetwork")

  # nodes carry one row per kept term; edges exist above the similarity cutoff
  expect_equal(nrow(w1$x$nodes), 5)
  expect_true(nrow(w1$x$edges) >= 1)
})

test_that("fig_enrich_visnet errors on empty input", {
  skip_if_not_installed("visNetwork")
  expect_error(fig_enrich_visnet(data.frame()), "No enrichment results")
})
