fake_gsea <- function() {
  data.frame(
    pathway = paste0("HALLMARK_", c("INFLAMMATORY_RESPONSE",
                     "INTERFERON_GAMMA_RESPONSE", "TNFA_SIGNALING_VIA_NFKB",
                     "E2F_TARGETS", "OXIDATIVE_PHOSPHORYLATION")),
    pval = 10^(-(20:16)), padj = 10^(-(18:14)),
    NES = c(2.9, 2.7, 2.5, -1.8, -1.2), ES = c(.7,.6,.6,-.5,-.4),
    size = c(200L, 190L, 195L, 180L, 200L),
    leadingEdge = I(list(paste0("g", 1:60), paste0("g", 30:90),
                         paste0("g", 1:55), paste0("g", 200:250),
                         paste0("g", 260:300))),
    stringsAsFactors = FALSE)
}

fake_ora <- function() {
  data.frame(
    ID = paste0("GO:", 1:4),
    Description = c("response to lipopolysaccharide", "response to virus",
                   "innate immune response", "cytokine production"),
    GeneRatio = c("60/400", "50/400", "45/400", "40/400"),
    BgRatio = "300/20000", pvalue = 10^(-(30:27)), padj = 10^(-(28:25)),
    qvalue = 10^(-(28:25)), Count = c(60L, 50L, 45L, 40L),
    geneID = c(paste(paste0("g", 1:60), collapse = "/"),
               paste(paste0("g", 30:79), collapse = "/"),
               paste(paste0("g", 1:45), collapse = "/"),
               paste(paste0("g", 100:139), collapse = "/")),
    stringsAsFactors = FALSE)
}

test_that("palette helpers return valid colors and scales", {
  cols <- rnaflow_colors(6)
  expect_length(cols, 6)
  expect_true(all(grepl("^#[0-9A-Fa-f]{6}$", cols)))
  expect_length(rnaflow_colors(30), 30)              # interpolates beyond base
  expect_true(all(grepl("^#", omics_ramp("batlow", 16))))
  expect_s3_class(scale_fill_omics_div("vik", limits = c(-2, 2)), "ScaleContinuous")
})

test_that("fig_enrich_map builds for GSEA and ORA", {
  skip_if_not_installed("ggraph")
  skip_if_not_installed("tidygraph")
  expect_s3_class(fig_enrich_map(fake_gsea(), n = 5), "ggplot")
  expect_s3_class(fig_enrich_map(fake_ora(), n = 4), "ggplot")
  expect_error(fig_enrich_map(fake_gsea()[0, ]), "No enrichment results")
})

test_that("fig_gsea_ridge builds from ranks, gene sets and GSEA table", {
  skip_if_not_installed("ggridges")
  set.seed(4)
  n <- 400
  res <- data.frame(gene = paste0("g", seq_len(n)),
                    log2FoldChange = rnorm(n), padj = runif(n),
                    pvalue = runif(n), stat = rnorm(n))
  gsea <- fake_gsea()
  gene_sets <- stats::setNames(gsea$leadingEdge, gsea$pathway)
  expect_s3_class(fig_gsea_ridge(res, gene_sets, gsea, n = 4), "ggplot")
})

test_that("fig_volcano glow adds a layer and still returns a ggplot", {
  res <- data.frame(gene = paste0("g", 1:200), baseMean = 50,
                    log2FoldChange = rnorm(200, 0, 2), lfcSE = 0.2,
                    stat = rnorm(200), pvalue = runif(200),
                    padj = pmin(runif(200) * 0.4, 1))
  p_plain <- fig_volcano(res)
  p_glow  <- fig_volcano(res, glow = TRUE)
  expect_s3_class(p_glow, "ggplot")
  expect_gt(length(p_glow$layers), length(p_plain$layers))
})

test_that("fig_module_network builds from a WGCNA result", {
  skip_if_not_installed("WGCNA")
  skip_if_not_installed("ggraph")
  set.seed(11)
  f1 <- rnorm(30); f2 <- rnorm(30)
  block <- function(f, k) t(sapply(seq_len(k),
    function(i) f * runif(1, .8, 1.2) + rnorm(length(f), 0, .25)))
  m <- rbind(block(f1, 40), block(f2, 40), matrix(rnorm(70 * 30), nrow = 70))
  rownames(m) <- paste0("g", seq_len(nrow(m)))
  colnames(m) <- paste0("s", seq_len(30))
  wg <- suppressWarnings(run_wgcna(wgcna_datexpr(m, 150), power = 6,
                                   min_module_size = 20))
  mod <- setdiff(unique(wg$modules), "grey")[1]
  expect_s3_class(fig_module_network(wg, mod, n = 25), "ggplot")
})
