# Compose RNAflow's own figures into three dense manuscript plates (Figures 1-3),
# six panels each (balanced 2x3 grids), full descriptive titles. Shows most of
# the tool's output. Caches the heavy objects. Run from the package root:
#   Rscript paper/make_panels.R
suppressPackageStartupMessages({
  if (!suppressWarnings(require(RNAflow, quietly = TRUE))) devtools::load_all(".", quiet = TRUE)
  library(ggplot2); library(patchwork)
})
set.seed(1)
ext   <- function(f) system.file("extdata", f, package = "RNAflow")
cache <- "paper/.panel_cache.rds"
OI <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#7F7F7F")

## ---- compute once, cache ----------------------------------------------
if (file.exists(cache)) { D <- readRDS(cache); message("loaded cache") } else {
  message("computing (first run, several minutes)...")
  ac <- read_counts(ext("demo_airway_counts.csv"))
  am <- read_metadata(ext("demo_airway_metadata.csv"), counts_samples = colnames(ac))
  a_res <- run_deseq2(ac, am, design = ~ cell + condition,
                      contrast = c("condition","Dex","Control"), shrink = TRUE)
  a_vst <- normalize_counts(ac, am, "vst")
  sets  <- get_gene_sets("human", collection = "H")
  a_gsea <- run_gsea(a_res, sets, rank_by = "stat")
  a_ora  <- run_ora(contrast_sig_genes(a_res, 0.05, 1), "human", db = "GO", ont = "BP",
                    universe = a_res$gene)
  tc <- read_counts(ext("demo_tcga_counts.csv"))
  tm <- read_metadata(ext("demo_tcga_metadata.csv"), counts_samples = colnames(tc))
  vst <- normalize_counts(tc, tm, "vst")
  pca <- compute_pca(vst, 500)
  gv  <- run_gsva(vst, sets, method = "gsva")
  datExpr <- wgcna_datexpr(vst, n_genes = 3500)
  sft <- wgcna_pick_power(datExpr); wg <- run_wgcna(datExpr, power = sft$suggested)
  mt  <- module_trait_cor(wg$MEs, build_traits(tm, rownames(datExpr)))
  mod_enrich <- tryCatch(enrich_modules(wg, "human", db = "GO", ont = "BP", n_per = 3),
                         error = function(e) NULL)
  dfs <- list()
  for (pr in list(c("LGG","LUAD"),c("KIRC","COAD"),c("BRCA","LUAD"),c("PRAD","THCA")))
    dfs[[paste(pr, collapse=" vs ")]] <- run_deseq2(tc, tm, design = ~ cancer_type,
      contrast = c("cancer_type", pr[1], pr[2]), shrink = FALSE)
  setlist <- contrast_sig_sets(dfs, 0.05, 1, "either")
  D <- list(am=am, a_res=a_res, a_vst=a_vst, a_gsea=a_gsea, a_ora=a_ora,
            tm=tm, vst=vst, pca=pca, gv=gv, sft=sft, wg=wg, mt=mt,
            mod_enrich=mod_enrich, dfs=dfs, setlist=setlist)
  saveRDS(D, cache)
}

## ---- theme + helpers --------------------------------------------------
theme_2026 <- function(b = 11.5) theme_publication() + theme(
  text = element_text(family = "sans"),
  plot.title = element_text(size = b - 0.5, face = "bold", margin = margin(b = 4)),
  axis.title.x = element_text(size = b - 2, margin = margin(t = 2)),
  axis.title.y = element_text(size = b - 2, margin = margin(r = 2)),
  axis.text = element_text(size = b - 3.5),
  legend.title = element_text(size = b - 3, face = "bold"),
  legend.text = element_text(size = b - 4), plot.margin = margin(7, 9, 7, 7))
tag_theme <- theme(plot.tag = element_text(size = 18, face = "bold", colour = "#1D9E75"))
gg <- function(p, title) p + labs(title = title) + theme_2026()
titled <- function(g, t) g + labs(title = t) +
  theme(plot.title = element_text(size = 11, face = "bold", margin = margin(b = 4)))
hm <- function(obj) {
  if (inherits(obj, "pheatmap")) ggplotify::as.ggplot(obj$gtable)
  else ggplotify::as.ggplot(grid::grid.grabExpr(
    if (inherits(obj, c("Heatmap","HeatmapList"))) ComplexHeatmap::draw(obj) else grid::grid.draw(obj)))
}
save_plate <- function(pl, name, w, h) {
  ggsave(paste0("paper/figures/", name, ".png"), pl, width = w, height = h, dpi = 300, bg = "white")
  ggsave(paste0("paper/figures/", name, ".pdf"), pl, width = w, height = h, bg = "white")
  cat("wrote", name, "\n")
}

## ---- Figure 1: differential expression & enrichment (airway) ----------
p1 <- wrap_plots(
  gg(fig_volcano(D$a_res, n_label = 8, mode = "publication"),
     "Differential expression: dexamethasone vs control (airway)"),
  gg(fig_ma(D$a_res, mode = "publication"), "MA plot of differential expression (airway)"),
  gg(fig_pval_hist(D$a_res, mode = "publication"), "P-value distribution (airway)"),
  titled(hm(fig_heatmap(D$a_vst, D$a_res, D$am, n_genes = 30, show_colnames = FALSE,
                        direction_annotation = TRUE, show_title = FALSE)),
         "Expression of top differential genes (airway)"),
  gg(fig_enrich_dot(D$a_gsea, n = 8), "Gene-set enrichment (airway; MSigDB Hallmark)"),
  gg(fig_enrich_bar(D$a_ora, n = 10, mode = "publication"),
     "Over-representation analysis (airway; GO BP)"),
  ncol = 3) + plot_annotation(tag_levels = "A") & tag_theme
save_plate(p1, "figure1", 17, 10.5)

## ---- Figure 2: structure & co-expression (TCGA) -----------------------
sc <- merge(D$pca$scores, D$tm, by = "sample")
p_pca <- ggplot(sc, aes(PC1, PC2, colour = cancer_type)) +
  geom_point(size = 2.4, alpha = .92, stroke = 0) +
  scale_colour_manual(values = OI, name = "Cancer type") +
  labs(title = "Principal-component analysis (TCGA pan-cancer cohort)",
       x = sprintf("PC1 (%.1f%%)", D$pca$pct[1]), y = sprintf("PC2 (%.1f%%)", D$pca$pct[2])) +
  theme_2026()
p2D <- if (!is.null(D$mod_enrich))
  (gg(fig_module_enrichment(D$mod_enrich, mode = "publication"),
      "Functional enrichment of co-expression modules (GO BP)") +
     theme(axis.text.x = ggtext::element_markdown(size = 7))) else
  gg(fig_module_sizes(D$wg, mode = "publication"), "WGCNA co-expression module sizes")
p2 <- wrap_plots(
  p_pca,
  titled(hm(fig_sample_cor(D$vst, D$tm, show_names = FALSE, title = NA)),
         "Sample-to-sample correlation (TCGA)"),
  titled(hm(fig_heatmap(D$vst, D$dfs[["LGG vs LUAD"]], D$tm, n_genes = 40,
                        show_colnames = FALSE, show_title = FALSE)),
         "Top differentially expressed genes across tumours (TCGA)"),
  gg(fig_soft_threshold(D$sft, mode = "publication"), "WGCNA scale-free soft-threshold selection"),
  gg(fig_module_trait(D$mt, mode = "publication"), "WGCNA module–trait correlation (TCGA)") +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 6.5)),
  p2D, ncol = 3) + plot_annotation(tag_levels = "A") & tag_theme
save_plate(p2, "figure2", 17, 11)

## ---- Figure 3: multi-contrast comparison & signatures (TCGA) ----------
p3 <- wrap_plots(
  gg(fig_volcano_grid(D$dfs, n_label = 3, mode = "publication"),
     "Pairwise differential expression across cancer types (TCGA)"),
  titled(hm(fig_upset(D$setlist, min_size = 1)),
         "Overlap of significant genes across contrasts (UpSet)"),
  titled(hm(fig_venn(D$setlist[1:3])),
         "Overlap of significant genes across contrasts (Venn)"),
  titled(hm(fig_lfc_heatmap(D$dfs, n_genes = 40, show_title = FALSE)),
         "Log2 fold-change of genes across contrasts"),
  gg(fig_contrast_alluvial(D$dfs, mode = "publication"),
     "Direction of change of genes across contrasts"),
  titled(hm(fig_gsva_heatmap(D$gv, D$tm, group_by = "cancer_type", n_top = 40, title = NA)),
         "Per-sample gene-set signatures (GSVA; Hallmark)"),
  ncol = 3) + plot_annotation(tag_levels = "A") & tag_theme
save_plate(p3, "figure3", 17, 11)
