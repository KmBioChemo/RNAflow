# Export every panel of Figures 2-4 as a standalone vector PDF (+ PNG) so the
# panels can be re-assembled by hand (e.g. in PowerPoint). Uses the cache.
# Run from the package root:  Rscript paper/make_panels_individual.R
suppressPackageStartupMessages({
  if (!suppressWarnings(require(RNAflow, quietly = TRUE))) devtools::load_all(".", quiet = TRUE)
  library(ggplot2)
})
set.seed(1)
cache <- "paper/.panel_cache.rds"
stopifnot(file.exists(cache)); D <- readRDS(cache)
sets  <- get_gene_sets("human", collection = "H")
outdir <- "paper/figures/panels"; dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
OI <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#7F7F7F")

theme_2026 <- function(b = 12) theme_publication() + theme(
  text = element_text(family = "sans"),
  plot.title = element_text(size = b, face = "bold", margin = margin(b = 5)),
  axis.title.x = element_text(size = b - 2, margin = margin(t = 2)),
  axis.title.y = element_text(size = b - 2, margin = margin(r = 2)),
  axis.text = element_text(size = b - 3), plot.margin = margin(9, 11, 9, 9))
gg <- function(p, title) p + labs(title = title) + theme_2026()

# save a ggplot / pheatmap / ComplexHeatmap / eulergram as its own PDF + PNG
savep <- function(obj, name, w, h, title = NULL) {
  base <- file.path(outdir, name)
  drawit <- function() {
    if (inherits(obj, "pheatmap")) { grid::grid.newpage(); grid::grid.draw(obj$gtable) }
    else if (inherits(obj, c("Heatmap", "HeatmapList"))) ComplexHeatmap::draw(obj)
    else { grid::grid.newpage(); grid::grid.draw(obj) }
  }
  if (inherits(obj, "ggplot")) {
    if (!is.null(title)) obj <- gg(obj, title)
    ggsave(paste0(base, ".pdf"), obj, width = w, height = h)
    ggsave(paste0(base, ".png"), obj, width = w, height = h, dpi = 300, bg = "white")
  } else {
    grDevices::pdf(paste0(base, ".pdf"), width = w, height = h); drawit(); grDevices::dev.off()
    grDevices::png(paste0(base, ".png"), width = w * 300, height = h * 300, res = 300); drawit(); grDevices::dev.off()
  }
  cat("wrote", name, "\n")
}

## ---- Figure 2 panels (airway: DE & enrichment) -----------------------
savep(gg(fig_volcano(D$a_res, n_label = 10, mode = "publication"),
         "Differential expression: dexamethasone vs control (airway)"), "F2A_volcano_airway", 6, 5)
savep(gg(fig_ma(D$a_res, mode = "publication"), "MA plot (airway)"), "F2B_ma_airway", 6, 5)
savep(gg(fig_pval_hist(D$a_res, mode = "publication"), "P-value distribution (airway)"), "F2C_pvalue_airway", 5.5, 4.5)
savep(fig_heatmap(D$a_vst, D$a_res, D$am, n_genes = 30, show_colnames = FALSE,
                  direction_annotation = TRUE, title = "Top differential genes (airway)"),
      "F2D_de_heatmap_airway", 6, 7)
savep(gg(fig_gsea_ridge(D$a_res, sets, D$a_gsea, n = 12, mode = "publication"),
         "Gene-set enrichment ridgeline (airway; MSigDB Hallmark)"), "F2E_gsea_ridge_airway", 7, 6)
savep(gg(fig_enrich_bar(D$a_ora, n = 12, mode = "publication"),
         "Over-representation analysis (airway; GO BP)"), "F2F_ora_bar_airway", 7.5, 5.5)

## ---- Figure 3 panels (TCGA: landscape & co-expression) ---------------
sc <- merge(D$pca$scores, D$tm, by = "sample")
p_pca <- ggplot(sc, aes(PC1, PC2, colour = cancer_type)) +
  geom_point(size = 2.8, alpha = .92, stroke = 0) +
  scale_colour_manual(values = OI, name = "Cancer type") +
  labs(title = "Principal-component analysis (TCGA)",
       x = sprintf("PC1 (%.1f%%)", D$pca$pct[1]), y = sprintf("PC2 (%.1f%%)", D$pca$pct[2])) +
  theme_2026()
savep(p_pca, "F3B_pca_tcga", 6.5, 5)
savep(fig_gsva_heatmap(D$gv, D$tm, group_by = "cancer_type", n_top = 45,
                       title = "Per-sample GSVA Hallmark signatures (TCGA)"),
      "F3A_gsva_heatmap_tcga", 9, 8)
savep(gg(fig_soft_threshold(D$sft, mode = "publication"),
         "WGCNA scale-free soft-threshold selection"), "F3C_wgcna_soft_threshold", 8, 4.8)
savep(gg(fig_module_trait(D$mt, mode = "publication"), "WGCNA module–trait correlation (TCGA)") +
        theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 8)),
      "F3D_wgcna_module_trait", 8, 6)
if (!is.null(D$mod_enrich))
  savep(gg(fig_module_enrichment(D$mod_enrich, max_terms = 20, mode = "publication"),
           "Functional enrichment of co-expression modules (GO BP)") +
          theme(axis.text.x = ggtext::element_markdown(size = 8)),
        "F3E_module_enrichment", 8.5, 7)

## ---- Figure 4 panels (TCGA: multi-contrast) --------------------------
savep(gg(fig_volcano_grid(D$dfs, n_label = 3, mode = "publication"),
         "Pairwise differential expression across cancer types (TCGA)"), "F4A_volcano_grid", 9, 7)
savep(fig_upset(D$setlist, min_size = 1), "F4B_upset", 8, 5.5)
savep(fig_venn(D$setlist[1:3]), "F4C_venn", 6.5, 6)
savep(fig_lfc_heatmap(D$dfs, n_genes = 40, title = "Log2 fold-change across contrasts"),
      "F4D_lfc_heatmap", 6.5, 8)
savep(gg(fig_contrast_alluvial(D$dfs, mode = "publication"),
         "Direction of change across contrasts"), "F4E_direction_alluvial", 7.5, 5)

cat("\nAll panels written to ", outdir, "\n")
