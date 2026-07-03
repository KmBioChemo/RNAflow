make_de <- function(genes, lfc, padj) {
  data.frame(
    gene = genes,
    baseMean = 100,
    log2FoldChange = lfc,
    lfcSE = 0.2,
    stat = lfc / 0.2,
    pvalue = padj,
    padj = padj,
    stringsAsFactors = FALSE
  )
}

# Two overlapping contrasts with hand-picked significance.
make_contrasts <- function() {
  list(
    A = make_de(c("g1", "g2", "g3", "g4"),
                lfc  = c( 2, -2,  0.1, 3),
                padj = c(0.001, 0.001, 0.4, 0.001)),   # sig: g1(up) g2(down) g4(up)
    B = make_de(c("g1", "g2", "g3", "g5"),
                lfc  = c( 2,  0.2, -3, -2),
                padj = c(0.001, 0.4, 0.001, 0.001))     # sig: g1(up) g3(down) g5(down)
  )
}

test_that("contrast_sig_genes respects thresholds and direction", {
  a <- make_contrasts()$A
  expect_setequal(contrast_sig_genes(a, 0.05, 1, "either"), c("g1", "g2", "g4"))
  expect_setequal(contrast_sig_genes(a, 0.05, 1, "up"),     c("g1", "g4"))
  expect_setequal(contrast_sig_genes(a, 0.05, 1, "down"),   c("g2"))
})

test_that("contrast_sig_genes drops NA padj / lfc", {
  d <- make_de(c("g1", "g2"), lfc = c(3, NA), padj = c(NA, 0.001))
  expect_length(contrast_sig_genes(d, 0.05, 1, "either"), 0)
})

test_that("contrast_sig_sets returns one named set per contrast", {
  sets <- contrast_sig_sets(make_contrasts(), 0.05, 1, "either")
  expect_named(sets, c("A", "B"))
  expect_setequal(sets$A, c("g1", "g2", "g4"))
  expect_setequal(sets$B, c("g1", "g3", "g5"))
})

test_that("contrast_lfc_matrix unions genes and fills missing with NA", {
  m <- contrast_lfc_matrix(make_contrasts())
  expect_equal(ncol(m), 2)
  expect_setequal(rownames(m), c("g1", "g2", "g3", "g4", "g5"))
  expect_equal(unname(m["g4", "A"]), 3)       # present in A
  expect_true(is.na(m["g4", "B"]))            # absent from B
  expect_true(is.na(m["g5", "A"]))            # absent from A
})

test_that("contrast_lfc_matrix honors a custom gene order", {
  m <- contrast_lfc_matrix(make_contrasts(), genes = c("g3", "g1"))
  expect_equal(rownames(m), c("g3", "g1"))
})

test_that("top_variable_genes keeps the most variable rows", {
  m <- matrix(c(0, 0,  5, -5,  1, 1.1),
              nrow = 3, byrow = TRUE,
              dimnames = list(c("flat", "swing", "mild"), c("A", "B")))
  out <- top_variable_genes(m, 2)
  expect_equal(nrow(out), 2)
  expect_true("swing" %in% rownames(out))
  expect_false("flat" %in% rownames(out))
})

test_that("check_contrasts rejects unnamed / duplicate / invalid input", {
  good <- make_contrasts()
  expect_error(contrast_sig_sets(unname(good)), "unique, non-empty name")
  expect_error(contrast_sig_sets(setNames(good, c("A", "A"))), "unique")
  bad <- list(X = data.frame(gene = "g1"))   # missing required columns
  expect_error(contrast_sig_sets(bad), "Contrast 'X'")
})

test_that("fig_volcano_grid returns a faceted ggplot", {
  p <- fig_volcano_grid(make_contrasts(), lfc_thr = 1, padj_thr = 0.05,
                        n_label = 2)
  expect_s3_class(p, "ggplot")
})

test_that("fig_lfc_heatmap returns a pheatmap object", {
  skip_if_not_installed("pheatmap")
  # widen the gene pool so >= 2 sig genes exist
  ct <- make_contrasts()
  hm <- fig_lfc_heatmap(ct, gene_src = "sig_union", padj_thr = 0.05, lfc_thr = 1)
  expect_s3_class(hm, "pheatmap")
})

test_that("fig_venn builds for 2-4 sets and rejects more", {
  skip_if_not_installed("eulerr")
  sets <- contrast_sig_sets(make_contrasts(), 0.05, 1, "either")
  expect_no_error(fig_venn(sets))
  big <- setNames(replicate(5, "g1", simplify = FALSE),
                  paste0("S", 1:5))
  expect_error(fig_venn(big), "fig_upset")
})

test_that("fig_upset builds a combination matrix plot", {
  skip_if_not_installed("ComplexHeatmap")
  sets <- contrast_sig_sets(make_contrasts(), 0.05, 1, "either")
  expect_no_error(fig_upset(sets))
})
