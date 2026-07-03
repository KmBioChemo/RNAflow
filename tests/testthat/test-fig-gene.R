make_expr <- function(n_genes = 50, n_samp = 12, seed = 3) {
  set.seed(seed)
  m <- matrix(rnorm(n_genes * n_samp, 8, 1.5), nrow = n_genes,
              dimnames = list(paste0("g", seq_len(n_genes)),
                              paste0("s", seq_len(n_samp))))
  md <- data.frame(sample = colnames(m),
                   group = rep(c("A", "B", "C"), length.out = n_samp),
                   stringsAsFactors = FALSE)
  list(m = m, md = md)
}

test_that("fig_gene_expression builds a ggplot for each style", {
  d <- make_expr()
  for (st in c("raincloud", "beeswarm", "box")) {
    p <- fig_gene_expression(d$m, d$md, gene = "g1", style = st)
    expect_s3_class(p, "ggplot")
  }
})

test_that("fig_gene_expression validates its inputs", {
  d <- make_expr()
  expect_error(fig_gene_expression(d$m, d$md, gene = "nope"), "not found")
  expect_error(
    fig_gene_expression(d$m, d$md[, 1, drop = FALSE], gene = "g1"),
    "annotation column")
})

mk_de <- function(genes, up, down, seed = 1) {
  set.seed(seed)
  n <- length(genes)
  lfc <- rnorm(n, 0, 0.3); padj <- runif(n, 0.2, 1)
  lfc[up] <- 2;  padj[up] <- 1e-4
  lfc[down] <- -2; padj[down] <- 1e-4
  data.frame(gene = genes, log2FoldChange = lfc, padj = padj,
             stringsAsFactors = FALSE)
}

test_that("contrast_direction_table keeps only genes significant somewhere", {
  genes <- paste0("g", 1:10)
  cs <- list(
    "c1" = mk_de(genes, up = 1:2, down = 3),
    "c2" = mk_de(genes, up = 2,   down = 8, seed = 2))
  tab <- contrast_direction_table(cs)
  expect_true(all(c("gene", "c1", "c2") %in% colnames(tab)))
  expect_true(all(levels(tab$c1) == c("Up", "NS", "Down")))
  # every kept gene is non-NS in at least one contrast
  nonNS <- apply(tab[, -1, drop = FALSE], 1, function(r) any(r != "NS"))
  expect_true(all(nonNS))
  expect_true("g1" %in% tab$gene)          # up in c1
  expect_false("g5" %in% tab$gene)         # NS everywhere
})

test_that("fig_contrast_alluvial builds a ggplot", {
  skip_if_not_installed("ggalluvial")
  genes <- paste0("g", 1:12)
  cs <- list("c1" = mk_de(genes, up = 1:3, down = 4:5),
             "c2" = mk_de(genes, up = 2:4, down = 6, seed = 5))
  expect_s3_class(fig_contrast_alluvial(cs), "ggplot")
})

test_that("contrast_direction_table needs >= 2 contrasts", {
  expect_error(contrast_direction_table(list(a = mk_de(paste0("g", 1:5), 1, 2))),
               "at least two")
})
