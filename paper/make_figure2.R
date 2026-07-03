# Figure 2 + downstream statistics for the RNAflow manuscript: the deeper
# TCGA analyses (co-expression modules, per-sample signatures, all-pairwise
# differential expression). Run from the package root:
#   Rscript paper/make_figure2.R
suppressPackageStartupMessages({
  if (!suppressWarnings(require(RNAflow, quietly = TRUE))) {
    devtools::load_all(".", quiet = TRUE)
  }
  library(ggplot2); library(patchwork)
})
set.seed(1)
ext <- function(f) system.file("extdata", f, package = "RNAflow")

t_counts <- read_counts(ext("demo_tcga_counts.csv"))
t_meta   <- read_metadata(ext("demo_tcga_metadata.csv"),
                          counts_samples = colnames(t_counts))
vst <- normalize_counts(t_counts, t_meta, "vst")
types <- sort(unique(t_meta$cancer_type))

tile_theme <- theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

## ---- Panel A: GSVA per-sample signatures, summarised by cancer type -----
pA <- tryCatch({
  sets <- get_gene_sets("human", collection = "H")
  gv <- run_gsva(vst, sets, method = "gsva")            # sets x samples
  # mean score per cancer type
  ct <- t_meta$cancer_type[match(colnames(gv), t_meta$sample)]
  m  <- sapply(types, function(k) rowMeans(gv[, ct == k, drop = FALSE]))
  # z-score each signature across types, keep the 18 most variable
  z  <- t(scale(t(m)))
  v  <- apply(z, 1, var); keep <- names(sort(v, decreasing = TRUE))[1:18]
  z  <- z[keep, ]
  cat(sprintf("[GSVA] %d Hallmark signatures scored across %d samples\n",
              nrow(gv), ncol(gv)))
  df <- data.frame(set = rep(rownames(z), ncol(z)),
                   type = rep(colnames(z), each = nrow(z)),
                   score = as.vector(z))
  df$set <- factor(df$set, levels = rev(keep))
  df$set <- factor(df$set, levels = rev(keep),
                   labels = rev(sub("^HALLMARK_", "", keep)))
  ggplot(df, aes(type, set, fill = score)) + geom_tile() +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                         name = "mean z") +
    labs(x = NULL, y = NULL, title = "TCGA: GSVA Hallmark signatures") +
    tile_theme
}, error = function(e) { cat("[GSVA panel skipped]:", conditionMessage(e), "\n"); NULL })

## ---- Panel B: WGCNA module-trait correlation ----------------------------
pB <- tryCatch({
  datExpr <- wgcna_datexpr(vst, n_genes = 3500)
  pw <- wgcna_pick_power(datExpr)$suggested
  wg <- run_wgcna(datExpr, power = pw)
  nmod <- length(setdiff(unique(wg$modules), "grey"))
  cat(sprintf("[WGCNA] soft-power=%s; %d modules (excl. grey); %d genes\n",
              pw, nmod, ncol(datExpr)))
  traits <- build_traits(t_meta, rownames(datExpr))
  mtc <- module_trait_cor(wg$MEs, traits)
  cormat <- if (is.list(mtc) && !is.null(mtc$cor)) mtc$cor else as.matrix(mtc)
  colnames(cormat) <- sub("^cancer_type[:_ ]*", "", colnames(cormat))
  rn <- sub("^ME", "", rownames(cormat))
  df <- data.frame(module = rep(rn, ncol(cormat)),
                   trait = rep(colnames(cormat), each = nrow(cormat)),
                   cor = as.vector(cormat))
  best <- df[which.max(abs(df$cor)), ]
  cat(sprintf("[WGCNA] strongest module-trait: %s ~ %s (r=%.2f)\n",
              best$module, best$trait, best$cor))
  ggplot(df, aes(trait, module, fill = cor)) + geom_tile() +
    scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                         limits = c(-1, 1), name = "r") +
    labs(x = NULL, y = "Module eigengene",
         title = "TCGA: WGCNA module-trait") + tile_theme
}, error = function(e) { cat("[WGCNA panel skipped]:", conditionMessage(e), "\n"); NULL })

## ---- Panel C: all-pairwise DE-gene counts -------------------------------
pC <- tryCatch({
  lst <- run_deseq2_all_pairs(t_counts, t_meta, design = ~ cancer_type,
                              shrink = FALSE)
  cat(sprintf("[all-pairwise] %d contrasts computed from one fit\n", length(lst)))
  cnt <- data.frame(a = character(), b = character(), n = integer())
  for (lab in names(lst)) {
    pr <- attr(lst[[lab]], "pair"); df <- lst[[lab]]
    nsig <- sum(!is.na(df$padj) & df$padj < 0.05 & abs(df$log2FoldChange) > 1)
    cnt <- rbind(cnt, data.frame(a = pr["treated"], b = pr["reference"], n = nsig),
                      data.frame(a = pr["reference"], b = pr["treated"], n = nsig))
  }
  cat(sprintf("[all-pairwise] sig genes per contrast: min=%d median=%d max=%d\n",
              min(cnt$n), as.integer(stats::median(cnt$n)), max(cnt$n)))
  cnt$a <- factor(cnt$a, levels = types); cnt$b <- factor(cnt$b, levels = types)
  ggplot(cnt, aes(a, b, fill = n)) + geom_tile() +
    scale_fill_viridis_c(option = "mako", direction = -1,
                         name = "DE genes") +
    labs(x = NULL, y = NULL, title = "TCGA: pairwise DE genes") + tile_theme
}, error = function(e) { cat("[pairwise panel skipped]:", conditionMessage(e), "\n"); NULL })

panels <- Filter(Negate(is.null), list(pA, pB, pC))
fig <- patchwork::wrap_plots(panels, nrow = 1) +
  patchwork::plot_annotation(tag_levels = "A")
ggsave("paper/figures/figure2.png", fig, width = 15, height = 5.2, dpi = 300,
       bg = "white")
ggsave("paper/figures/figure2.pdf", fig, width = 15, height = 5.2, bg = "white")
cat("Wrote paper/figures/figure2.{png,pdf}\n")
