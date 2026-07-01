test_that("run_deseq2 keeps the Wald stat and separates inference from shrinkage", {
  skip_if_not_installed("DESeq2")
  skip_on_cran()
  f_counts <- system.file("extdata", "demo_counts.csv", package = "RNAflow")
  f_meta   <- system.file("extdata", "demo_metadata.csv", package = "RNAflow")
  skip_if(!file.exists(f_counts) || !file.exists(f_meta))
  counts <- read_counts(f_counts)
  meta   <- read_metadata(f_meta, counts_samples = colnames(counts))

  ct <- c("condition", "Treatment", "Control")
  raw <- suppressMessages(suppressWarnings(
    run_deseq2(counts, meta, design = ~condition, contrast = ct, shrink = FALSE)))

  # apeglm actually engages here (Control is the reference level)
  shr_type <- if (requireNamespace("apeglm", quietly = TRUE)) "apeglm" else "normal"
  shk <- suppressMessages(suppressWarnings(
    run_deseq2(counts, meta, design = ~condition, contrast = ct,
               shrink = TRUE, shrink_type = shr_type)))

  # C1: a `stat` column is present even under apeglm (which normally drops it)
  expect_true("stat" %in% colnames(shk))
  expect_false(any(is.na(match(c("stat", "pvalue", "padj"), colnames(shk)))))

  # Inference (stat / pvalue / padj) is identical to the unshrunken test
  expect_equal(shk$stat, raw$stat)
  expect_equal(shk$padj, raw$padj)

  # Effect size is shrunken (magnitudes no larger than the MLE, on average)
  expect_lte(mean(abs(shk$log2FoldChange), na.rm = TRUE),
             mean(abs(raw$log2FoldChange), na.rm = TRUE) + 1e-8)

  # The estimator actually used is recorded
  expect_true(attr(shk, "shrink") %in% c("apeglm", "normal"))

  # And GSEA ranking by the default metric now works on shrunken results
  expect_gt(length(rank_genes(shk, by = "stat")), 0)
})

test_that("run_deseq2 keeps numeric covariates numeric (continuous adjustment)", {
  skip_if_not_installed("DESeq2")
  skip_on_cran()
  f_counts <- system.file("extdata", "demo_counts.csv", package = "RNAflow")
  f_meta   <- system.file("extdata", "demo_metadata.csv", package = "RNAflow")
  skip_if(!file.exists(f_counts) || !file.exists(f_meta))
  counts <- read_counts(f_counts)
  meta   <- read_metadata(f_meta, counts_samples = colnames(counts))

  # A continuous covariate, distinct per sample. If it were coerced to a
  # factor (one level per sample) the design would be rank-deficient and
  # DESeq2 would error; keeping it numeric makes it a valid adjustment.
  meta$rin <- seq(4, 9, length.out = nrow(meta))

  res <- suppressMessages(suppressWarnings(
    run_deseq2(counts, meta, design = ~ rin + condition,
               contrast = c("condition", "Treatment", "Control"),
               shrink = FALSE)))
  expect_s3_class(res, "data.frame")
  expect_gt(nrow(res), 0)
  expect_true(all(c("stat", "padj") %in% colnames(res)))
})
