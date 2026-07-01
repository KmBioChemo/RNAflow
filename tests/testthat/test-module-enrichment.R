fake_combined <- function() {
  do.call(rbind, lapply(c("blue", "brown", "yellow"), function(m) {
    data.frame(
      ID = paste0("GO:", m, 1:4),
      Description = paste(m, "term", 1:4),
      GeneRatio = "20/300", BgRatio = "100/20000",
      pvalue = 10^(-(10:7)), padj = 10^(-(9:6)), qvalue = 10^(-(9:6)),
      Count = c(40L, 30L, 25L, 20L), geneID = "a/b/c",
      module = m, stringsAsFactors = FALSE)
  }))
}

test_that("fig_module_enrichment returns a ggplot", {
  expect_s3_class(fig_module_enrichment(fake_combined()), "ggplot")
  expect_s3_class(fig_module_enrichment(fake_combined(), mode = "publication"),
                  "ggplot")
})

test_that("fig_module_enrichment respects max_terms", {
  p <- fig_module_enrichment(fake_combined(), max_terms = 5)
  # data behind the plot keeps at most 5 distinct terms
  expect_lte(length(unique(p$data$Description)), 5)
})

test_that("fig_module_enrichment validates its input", {
  expect_error(fig_module_enrichment(data.frame(x = 1)), "enrich_modules")
  expect_error(fig_module_enrichment(fake_combined()[0, ]), "enrich_modules")
})

test_that("enrich_modules runs per module on a real WGCNA object", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Mm.eg.db")
  skip_on_cran()
  f_counts <- system.file("extdata", "demo_multi_counts.csv", package = "RNAflow")
  f_meta   <- system.file("extdata", "demo_multi_metadata.csv", package = "RNAflow")
  skip_if(!file.exists(f_counts) || !file.exists(f_meta))
  counts <- read_counts(f_counts)
  meta   <- read_metadata(f_meta, counts_samples = colnames(counts))
  norm <- suppressMessages(suppressWarnings(
    normalize_counts(counts, meta, "vst")))
  wg <- suppressWarnings(run_wgcna(wgcna_datexpr(norm, 1500), power = 12,
                                   min_module_size = 30))
  comb <- suppressMessages(suppressWarnings(
    enrich_modules(wg, "mouse", db = "GO", ont = "BP", n_per = 3)))
  expect_true("module" %in% colnames(comb))
  expect_gt(nrow(comb), 0)
  expect_s3_class(fig_module_enrichment(comb), "ggplot")
})
