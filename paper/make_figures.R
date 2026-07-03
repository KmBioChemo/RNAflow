# Generate the manuscript figures for RNAflow *using the package's own figure
# functions* (fig_volcano, fig_enrich_dot, fig_module_trait, fig_gsva_heatmap,
# fig_upset, ...), on the two bundled datasets, and print the exact statistics
# quoted in the Results. Figures are therefore the tool's real outputs.
# Run from the package root:  Rscript paper/make_figures.R
suppressPackageStartupMessages({
  if (!suppressWarnings(require(RNAflow, quietly = TRUE))) {
    devtools::load_all(".", quiet = TRUE)
  }
  library(ggplot2); library(patchwork)
})
set.seed(1)
OI  <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
         "#0072B2", "#D55E00", "#CC79A7", "#7F7F7F")
ext <- function(f) system.file("extdata", f, package = "RNAflow")
dir.create("paper/figures", showWarnings = FALSE, recursive = TRUE)

## ======================================================================
## airway: simple two-group study (clean DE + enrichment)
## ======================================================================
a_counts <- read_counts(ext("demo_airway_counts.csv"))
a_meta   <- read_metadata(ext("demo_airway_metadata.csv"),
                          counts_samples = colnames(a_counts))
a_res <- run_deseq2(a_counts, a_meta, design = ~ cell + condition,
                    contrast = c("condition", "Dex", "Control"), shrink = TRUE)
a_sig <- a_res[!is.na(a_res$padj) & a_res$padj < 0.05 & abs(a_res$log2FoldChange) > 1, ]
cat(sprintf("[airway] %d tested; %d sig (%d up, %d down)\n",
            nrow(a_res), nrow(a_sig),
            sum(a_sig$log2FoldChange > 0), sum(a_sig$log2FoldChange < 0)))

sets  <- get_gene_sets("human", collection = "H")
a_gsea <- run_gsea(a_res, sets, rank_by = "stat")
cat(sprintf("[airway] GSEA Hallmark: %d sets padj<0.05\n",
            sum(a_gsea$padj < 0.05, na.rm = TRUE)))

p_volcano <- fig_volcano(a_res, lfc_thr = 1, padj_thr = 0.05, n_label = 12,
                         mode = "publication") +
  ggtitle("airway: dexamethasone vs control")
p_gsea <- fig_enrich_dot(a_gsea, n = 10) + ggtitle("airway: GSEA (Hallmark)")

## ======================================================================
## TCGA: complex eight-cancer-type cohort (the "big" dataset)
## ======================================================================
t_counts <- read_counts(ext("demo_tcga_counts.csv"))
t_meta   <- read_metadata(ext("demo_tcga_metadata.csv"),
                          counts_samples = colnames(t_counts))
vst   <- normalize_counts(t_counts, t_meta, "vst")
types <- sort(unique(t_meta$cancer_type))

# -- PCA (the tool's PCA engine; rendered static for print) --
pca <- compute_pca(vst, 500)
sc  <- merge(pca$scores, t_meta, by = "sample")
cat(sprintf("[TCGA] PCA PC1=%.1f%% PC2=%.1f%%\n", pca$pct[1], pca$pct[2]))
p_pca <- ggplot(sc, aes(PC1, PC2, colour = cancer_type)) +
  geom_point(size = 2.4, alpha = 0.9) +
  scale_colour_manual(values = OI, name = "Cancer type") +
  labs(x = sprintf("PC1 (%.1f%%)", pca$pct[1]),
       y = sprintf("PC2 (%.1f%%)", pca$pct[2]), title = "TCGA: PCA") +
  theme_publication()

# -- WGCNA module-trait (the tool's fig_module_trait) --
datExpr <- wgcna_datexpr(vst, n_genes = 3500)
pw <- wgcna_pick_power(datExpr)$suggested
wg <- run_wgcna(datExpr, power = pw)
nmod <- length(setdiff(unique(wg$modules), "grey"))
traits <- build_traits(t_meta, rownames(datExpr))
mt <- module_trait_cor(wg$MEs, traits)
cat(sprintf("[TCGA] WGCNA power=%s, %d modules\n", pw, nmod))
p_mt <- fig_module_trait(mt, mode = "publication") +
  ggtitle("TCGA: WGCNA module-trait")

## ---- Figure 1: DE + exploration + enrichment + networks --------------
fig1 <- (p_pca | p_volcano) / (p_gsea | p_mt) +
  patchwork::plot_annotation(tag_levels = "A")
ggsave("paper/figures/figure1.png", fig1, width = 13, height = 10, dpi = 300, bg = "white")
ggsave("paper/figures/figure1.pdf", fig1, width = 13, height = 10, bg = "white")
cat("Wrote figure1\n")

## ---- Figure 2: per-sample GSVA signatures (the tool's fig_gsva_heatmap)
gv <- run_gsva(vst, sets, method = "gsva")
cat(sprintf("[TCGA] GSVA: %d signatures x %d samples\n", nrow(gv), ncol(gv)))
ph_gsva <- fig_gsva_heatmap(gv, t_meta, group_by = "cancer_type", n_top = 40,
                            title = "TCGA: per-sample GSVA Hallmark signatures")
grDevices::png("paper/figures/figure2.png", width = 3000, height = 2400, res = 300)
grid::grid.newpage(); grid::grid.draw(ph_gsva$gtable); grDevices::dev.off()
grDevices::pdf("paper/figures/figure2.pdf", width = 10, height = 8)
grid::grid.newpage(); grid::grid.draw(ph_gsva$gtable); grDevices::dev.off()
cat("Wrote figure2\n")

## ---- Figure 3: multi-contrast overlap (the tool's fig_upset) ----------
pairs <- list(c("LUAD","LGG"), c("BRCA","LUAD"), c("KIRC","COAD"),
              c("PRAD","THCA"), c("SKCM","BRCA"))
setlist <- list()
for (pr in pairs) {
  r <- run_deseq2(t_counts, t_meta, design = ~ cancer_type,
                  contrast = c("cancer_type", pr[1], pr[2]), shrink = FALSE)
  sg <- r$gene[!is.na(r$padj) & r$padj < 0.05 & abs(r$log2FoldChange) > 1]
  setlist[[paste(pr, collapse = " vs ")]] <- sg
}
cat(sprintf("[TCGA] UpSet over %d contrasts; sizes: %s\n", length(setlist),
            paste(vapply(setlist, length, 0L), collapse = ", ")))
up <- fig_upset(setlist, min_size = 1)
grDevices::png("paper/figures/figure3.png", width = 3200, height = 1900, res = 300)
ComplexHeatmap::draw(up); grDevices::dev.off()
grDevices::pdf("paper/figures/figure3.pdf", width = 11, height = 6.5)
ComplexHeatmap::draw(up); grDevices::dev.off()
cat("Wrote figure3\n")
