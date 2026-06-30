make_ranked_de <- function(n = 1000, seed = 7) {
  set.seed(seed)
  data.frame(
    gene = paste0("G", seq_len(n)),
    baseMean = 100,
    log2FoldChange = rnorm(n),
    lfcSE = 0.2,
    stat = rnorm(n),
    pvalue = runif(n),
    padj = runif(n),
    stringsAsFactors = FALSE
  )
}

test_that("rank_genes produces a sorted, de-duplicated named vector", {
  res <- make_ranked_de(50)
  r <- rank_genes(res, by = "stat")
  expect_type(r, "double")
  expect_false(is.unsorted(rev(r)))          # decreasing
  expect_equal(length(r), 50)
  expect_true(all(names(r) %in% res$gene))
})

test_that("rank_genes by signed_p uses p-value and fold-change sign", {
  res <- data.frame(gene = c("a", "b"), log2FoldChange = c(2, -2),
                    pvalue = c(1e-4, 1e-4), padj = c(1e-3, 1e-3),
                    stat = c(5, -5))
  r <- rank_genes(res, by = "signed_p")
  expect_gt(r["a"], 0)
  expect_lt(r["b"], 0)
})

test_that("run_gsea recovers a planted enriched set", {
  skip_if_not_installed("fgsea")
  res <- make_ranked_de(1000)
  res <- res[order(res$stat, decreasing = TRUE), ]
  planted <- head(res$gene, 60)               # the most up-ranked genes
  random  <- sample(res$gene, 60)
  gs <- list(PLANTED_UP = planted, RANDOM = random)

  out <- run_gsea(res, gs, rank_by = "stat", min_size = 15)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("pathway", "padj", "NES") %in% colnames(out)))
  planted_row <- out[out$pathway == "PLANTED_UP", ]
  expect_gt(planted_row$NES, 1.5)             # strongly positive
  expect_lt(planted_row$pval, 0.05)
})

test_that("run_gsea validates inputs", {
  res <- make_ranked_de(50)
  expect_error(run_gsea(res, list()), "non-empty named list")
})

# ---- integration with the real annotation stack -----------------------------

test_that("get_gene_sets returns Hallmark sets for mouse", {
  skip_if_not_installed("msigdbr")
  gs <- get_gene_sets("mouse", collection = "H")
  expect_type(gs, "list")
  expect_gt(length(gs), 40)                   # ~50 Hallmark sets
  expect_true(all(grepl("^HALLMARK_", names(gs))))
  expect_type(gs[[1]], "character")
})

test_that("run_ora on GO BP returns a tidy frame for mouse symbols", {
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Mm.eg.db")
  # A handful of real mouse symbols; result may be empty, but must be tidy.
  genes <- c("Actb", "Gapdh", "Tnf", "Il6", "Il1b", "Nfkb1", "Stat1",
             "Cxcl10", "Ccl2", "Irf7")
  out <- run_ora(genes, "mouse", db = "GO", ont = "BP")
  expect_s3_class(out, "data.frame")
  expect_true(all(c("ID", "Description", "padj", "Count") %in% colnames(out)))
})
