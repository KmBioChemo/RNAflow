test_that("run_deseq2_all_pairs returns every pairwise contrast from one fit", {
  skip_if_not_installed("DESeq2")
  skip_on_cran()
  set.seed(7)
  groups  <- rep(c("A", "B", "C"), each = 5)          # 15 samples, 3 groups
  n_genes <- 300
  base <- matrix(rnbinom(n_genes * length(groups), mu = 100, size = 10),
                 nrow = n_genes)
  base[1:20,  groups == "B"] <- base[1:20,  groups == "B"] * 4   # B-up block
  base[21:40, groups == "C"] <- base[21:40, groups == "C"] * 4   # C-up block
  storage.mode(base) <- "integer"
  rownames(base) <- paste0("g", seq_len(n_genes))
  colnames(base) <- paste0("s", seq_along(groups))
  meta <- data.frame(sample = colnames(base), grp = groups,
                     stringsAsFactors = FALSE)

  res <- suppressMessages(suppressWarnings(
    run_deseq2_all_pairs(base, meta, design = ~ grp, shrink = FALSE)))

  expect_length(res, 3)                                # choose(3, 2)
  expect_setequal(names(res),
                  c("grp: B vs A", "grp: C vs A", "grp: C vs B"))
  for (df in res) {
    expect_s3_class(df, "data.frame")
    expect_true(all(c("gene", "log2FoldChange", "padj", "stat") %in% colnames(df)))
    pr <- attr(df, "pair")
    expect_true(all(c("treated", "reference") %in% names(pr)))
  }
  # the planted B-up block is significant and up in "B vs A"
  ba <- res[["grp: B vs A"]]
  expect_gte(sum(ba$padj < 0.05 & ba$log2FoldChange > 1, na.rm = TRUE), 5)
})

test_that("run_deseq2_all_pairs needs at least two levels", {
  skip_if_not_installed("DESeq2")
  base <- matrix(rpois(300, 50), nrow = 30,
                 dimnames = list(paste0("g", 1:30), paste0("s", 1:10)))
  storage.mode(base) <- "integer"
  meta <- data.frame(sample = colnames(base), grp = "A")
  expect_error(run_deseq2_all_pairs(base, meta, design = ~ grp),
               "fewer than 2 levels")
})
