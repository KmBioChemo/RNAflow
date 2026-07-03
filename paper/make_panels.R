# Compose RNAflow's own figures into the three imposing manuscript plates
# (Figures 1-3), following FIGURE_PLAN.md. Computes once and caches the heavy
# objects, so layout tweaks re-run fast. Run from the package root:
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
  message("computing (first run, a few minutes)...")
  ac <- read_counts(ext("demo_airway_counts.csv"))
  am <- read_metadata(ext("demo_airway_metadata.csv"), counts_samples = colnames(ac))
  a_res <- run_deseq2(ac, am, design = ~ cell + condition,
                      contrast = c("condition","Dex","Control"), shrink = TRUE)
  sets  <- get_gene_sets("human", collection = "H")
  a_gsea <- run_gsea(a_res, sets, rank_by = "stat")
  tc <- read_counts(ext("demo_tcga_counts.csv"))
  tm <- read_metadata(ext("demo_tcga_metadata.csv"), counts_samples = colnames(tc))
  vst <- normalize_counts(tc, tm, "vst")
  pca <- compute_pca(vst, 500)
  gv  <- run_gsva(vst, sets, method = "gsva")
  datExpr <- wgcna_datexpr(vst, n_genes = 3500)
  sft <- wgcna_pick_power(datExpr); wg <- run_wgcna(datExpr, power = sft$suggested)
  mt  <- module_trait_cor(wg$MEs, build_traits(tm, rownames(datExpr)))
  grp <- tm$cancer_type[match(rownames(datExpr), tm$sample)]
  dfs <- list()
  for (pr in list(c("LGG","LUAD"),c("KIRC","COAD"),c("BRCA","LUAD"),c("PRAD","THCA")))
    dfs[[paste(pr, collapse=" vs ")]] <- run_deseq2(tc, tm, design = ~ cancer_type,
      contrast = c("cancer_type", pr[1], pr[2]), shrink = FALSE)
  setlist <- contrast_sig_sets(dfs, 0.05, 1, "either")
  D <- list(a_res=a_res, a_gsea=a_gsea, tm=tm, pca=pca, gv=gv, sft=sft, wg=wg,
            mt=mt, grp=grp, dfs=dfs, setlist=setlist)
  saveRDS(D, cache)
}

## ---- "2026" theme + helpers -------------------------------------------
theme_2026 <- function(b = 12) theme_publication() + theme(
  text = element_text(family = "sans"),
  plot.title = element_text(size = b, face = "bold", margin = margin(b = 4)),
  plot.subtitle = element_text(size = b - 2, colour = "#5a6472"),
  axis.title = element_text(size = b - 2), axis.text = element_text(size = b - 3.5),
  legend.title = element_text(size = b - 3, face = "bold"),
  legend.text = element_text(size = b - 4), plot.margin = margin(10, 12, 10, 12))
tag_theme <- theme(plot.tag = element_text(size = 22, face = "bold", colour = "#1D9E75"))
gg <- function(p, title, sub = NULL) p + labs(title = title, subtitle = sub) + theme_2026()
hm <- function(obj) {                                    # heatmap/upset/venn -> ggplot
  if (inherits(obj, "pheatmap")) ggplotify::as.ggplot(obj$gtable)
  else ggplotify::as.ggplot(grid::grid.grabExpr(
    if (inherits(obj, c("Heatmap","HeatmapList"))) ComplexHeatmap::draw(obj) else grid::grid.draw(obj)))
}
save_plate <- function(pl, name, w, h) {
  ggsave(paste0("paper/figures/", name, ".png"), pl, width = w, height = h, dpi = 300, bg = "white")
  ggsave(paste0("paper/figures/", name, ".pdf"), pl, width = w, height = h, bg = "white")
  cat("wrote", name, "\n")
}

## ---- Plate 1 (Figure 1): overview, hero = PCA -------------------------
sc <- merge(D$pca$scores, D$tm, by = "sample")
p_pca <- ggplot(sc, aes(PC1, PC2, colour = cancer_type)) +
  geom_point(size = 3, alpha = .92, stroke = 0) +
  scale_colour_manual(values = OI, name = "Cancer type") +
  labs(title = "TCGA pan-cancer: PCA", subtitle = "120 tumours, 8 cancer types",
       x = sprintf("PC1 (%.1f%%)", D$pca$pct[1]), y = sprintf("PC2 (%.1f%%)", D$pca$pct[2])) +
  theme_2026(13)
p1 <- wrap_plots(
  gg(fig_volcano(D$a_res, n_label = 10, mode = "publication"), "airway: DE", "dexamethasone vs control"),
  gg(fig_enrich_dot(D$a_gsea, n = 8), "airway: GSEA", "MSigDB Hallmark"),
  p_pca,
  gg(fig_module_trait(D$mt, mode = "publication"), "TCGA: WGCNA module-trait"),
  design = "AC\nBC\nDD", heights = c(1, 1, 1.05)) +
  plot_annotation(tag_levels = "A") & tag_theme
save_plate(p1, "figure1", 16, 13)

## ---- Plate 2 (Figure 2): signatures + co-expression, hero = GSVA ------
p2 <- wrap_plots(
  hm(fig_gsva_heatmap(D$gv, D$tm, group_by = "cancer_type", n_top = 40,
                      title = "TCGA: per-sample GSVA Hallmark signatures")),
  gg(fig_soft_threshold(D$sft, mode = "publication"), "WGCNA soft-threshold"),
  gg(fig_module_sizes(D$wg, mode = "publication"), "WGCNA module sizes"),
  gg(fig_eigengene(D$wg, "turquoise", groups = D$grp, mode = "publication"),
     "Turquoise eigengene by cancer type") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          axis.title.x = element_blank()),
  design = "AAB\nAAC\nDDD", heights = c(1, 1, 0.9)) +
  plot_annotation(tag_levels = "A") & tag_theme
save_plate(p2, "figure2", 16, 12)

## ---- Plate 3 (Figure 3): multi-contrast, hero = volcano grid ----------
p3 <- wrap_plots(
  gg(fig_volcano_grid(D$dfs, n_label = 3, mode = "publication"), "TCGA: pairwise volcano grid"),
  hm(fig_upset(D$setlist, min_size = 1)),
  hm(fig_venn(D$setlist[1:3])),
  hm(fig_lfc_heatmap(D$dfs, n_genes = 40)),
  gg(fig_contrast_alluvial(D$dfs, mode = "publication"), "Direction alluvial"),
  design = "AAB\nAAC\nDDE", heights = c(1, 1, 0.95)) +
  plot_annotation(tag_levels = "A") & tag_theme
save_plate(p3, "figure3", 16, 12)
