# Compose RNAflow's own figure outputs into the three composed manuscript plates:
#   figure2 = differential expression & enrichment (airway, 2x3 grid)
#   figure3 = molecular characterisation (TCGA; GSVA hero + PCA + WGCNA)
#   figure4 = multi-contrast comparison (TCGA; volcano-grid hero + overlap/flow)
# See paper/FIGURE_PLAN.md for the panel->figure map and paper.md for the
# captions (this script must stay in sync with both). figure1 (overview) and
# figure5 (validation) are produced elsewhere (GRAPHICAL_ABSTRACT.md /
# make_validation.R). Caches the heavy objects. Run from the package root:
#   Rscript paper/make_panels.R
suppressPackageStartupMessages({
  if (!suppressWarnings(require(RNAflow, quietly = TRUE))) devtools::load_all(".", quiet = TRUE)
  library(ggplot2); library(patchwork); library(ragg); library(png)
})
set.seed(1)
# Portable Helvetica-metric font: Helvetica/Arial on macOS/Windows, Liberation Sans on Linux.
.pick_font <- function() {
  pref <- c("Helvetica", "Arial", "Liberation Sans", "DejaVu Sans")
  fams <- tryCatch(systemfonts::system_fonts()$family, error = function(e) character(0))
  hit <- pref[pref %in% fams]; if (length(hit)) hit[1] else "sans"
}
FONT <- .pick_font()
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
  text = element_text(family = FONT),
  plot.title = element_text(size = b - 0.5, face = "bold", margin = margin(b = 4)),
  axis.title.x = element_text(size = b - 2, margin = margin(t = 2)),
  axis.title.y = element_text(size = b - 2, margin = margin(r = 2)),
  axis.text = element_text(size = b - 3.5),
  legend.title = element_text(size = b - 3, face = "bold"),
  legend.text = element_text(size = b - 4), plot.margin = margin(7, 9, 7, 7))
tag_theme <- theme(plot.tag = element_text(size = 16, face = "bold", colour = "black"))
gg <- function(p, title) p + labs(title = title) + theme_2026()
titled <- function(g, t) g + labs(title = t) +
  theme(plot.title = element_text(size = 11, face = "bold", margin = margin(b = 4)))
# Render a heatmap object (pheatmap / ComplexHeatmap / grid grob) to a PNG sized
# to its target cell, then embed it FILLING the panel. pheatmap & ComplexHeatmap
# reflow to fill whatever device size they are drawn to, so this uses all the
# available space with no letterboxing and no deformation (unlike as.ggplot(),
# which preserves the object's own fixed-unit layout and leaves white margins).
hm_panel <- function(obj, w, h, title = NULL, dpi = 320) {
  f <- tempfile(fileext = ".png")
  ragg::agg_png(f, width = w, height = h, units = "in", res = dpi, background = "white")
  if (inherits(obj, "pheatmap")) grid::grid.draw(obj$gtable)
  else if (inherits(obj, c("Heatmap", "HeatmapList"))) ComplexHeatmap::draw(obj)
  else grid::grid.draw(obj)
  grDevices::dev.off()
  img <- png::readPNG(f)
  g <- ggplot() +
    annotation_custom(grid::rasterGrob(img, width = unit(1, "npc"), height = unit(1, "npc"))) +
    coord_cartesian(expand = FALSE) + theme_void() + theme(plot.margin = margin(2, 2, 2, 2))
  if (!is.null(title)) g <- g + ggtitle(title) +
    theme(plot.title = element_text(size = 11, face = "bold", family = FONT,
                                    hjust = 0, margin = margin(b = 3)))
  g
}
save_plate <- function(pl, name, w, h) {
  ggsave(paste0("paper/figures/", name, ".png"), pl, width = w, height = h, dpi = 300,
         bg = "white", device = ragg::agg_png)
  ggsave(paste0("paper/figures/", name, ".pdf"), pl, width = w, height = h, bg = "white",
         device = grDevices::cairo_pdf)
  cat("wrote", name, "\n")
}

## ---- figure2: differential expression & enrichment (airway) -----------
# grouped: DE diagnostics (volcano/MA/p-value) then enrichment (GSEA ridge/ORA)
sets <- get_gene_sets("human", collection = "H")
# 16 x 10, 2x3 grid -> each cell ~ 5.33w x 5.0h; heatmap rendered to fill its cell (minus title)
pD2 <- hm_panel(fig_heatmap(D$a_vst, D$a_res, D$am, n_genes = 30, show_colnames = FALSE,
                            direction_annotation = TRUE, show_title = FALSE),
                w = 5.3, h = 4.55, title = "Top differential genes")
p1 <- wrap_plots(
  gg(fig_volcano(D$a_res, n_label = 8, mode = "publication"), "Differential expression"),
  gg(fig_ma(D$a_res, mode = "publication"), "MA plot"),
  gg(fig_pval_hist(D$a_res, mode = "publication"), "P-value histogram"),
  pD2,
  gg(fig_gsea_ridge(D$a_res, sets, D$a_gsea, n = 10, mode = "publication"),
     "Gene-set enrichment (GSEA)"),
  gg(fig_enrich_bar(D$a_ora, n = 10, mode = "publication"),
     "Over-representation (GO BP)"),
  ncol = 3) + plot_annotation(tag_levels = "A") & tag_theme
save_plate(p1, "figure2", 16, 10)

## ---- figure3: molecular landscape & co-expression (TCGA) --------------
# hero: the per-sample GSVA signature heatmap; supporting: PCA + WGCNA
sc <- merge(D$pca$scores, D$tm, by = "sample")
p_pca <- ggplot(sc, aes(PC1, PC2, colour = cancer_type)) +
  geom_point(size = 2.4, alpha = .92, stroke = 0) +
  scale_colour_manual(values = OI, name = "Cancer type") +
  labs(title = "Principal-component analysis",
       x = sprintf("PC1 (%.1f%%)", D$pca$pct[1]), y = sprintf("PC2 (%.1f%%)", D$pca$pct[2])) +
  theme_2026()
p2E <- if (!is.null(D$mod_enrich))
  (gg(fig_module_enrichment(D$mod_enrich, max_terms = 15, mode = "publication"),
      "Module enrichment (GO BP)") +
     theme(axis.text.x = ggtext::element_markdown(size = 7))) else
  gg(fig_module_sizes(D$wg, mode = "publication"), "Module sizes")
# 16.5 x 9.5, design AABC/AADE (4 cols, 2 rows): A = GSVA hero (cols1-2, both rows)
# -> ~8.0w x 9.0h; rendered to fill that tall cell. B..E are single ~4.0 x 4.5 cells.
pA3 <- hm_panel(fig_gsva_heatmap(D$gv, D$tm, group_by = "cancer_type", n_top = 45, title = ""),
                w = 7.9, h = 8.6, title = "Per-sample GSVA signatures (Hallmark)")
p2 <- wrap_plots(
  pA3,
  p_pca,
  gg(fig_soft_threshold(D$sft, mode = "publication"), "Soft-threshold selection"),
  gg(fig_module_trait(D$mt, mode = "publication"), "Module-trait correlation") +
    theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 6.5)),
  p2E,
  # GSVA heatmap needs many rows -> tall hero on the left; four uniform-height panels on the right
  design = "AABC\nAADE") + plot_annotation(tag_levels = "A") & tag_theme
save_plate(p2, "figure3", 16.5, 9.5)

## ---- figure4: multi-contrast comparison (TCGA) ------------------------
# hero: the pairwise volcano grid; supporting: overlap + flow views
# 16 x 9.5, design AAABBB/CCDDEE (6 cols, 2 rows): top row A|B halves (~8 x 4.6),
# bottom row C|D|E thirds (~5.3 x 4.6). Overlap heatmaps rendered to fill their cells.
pB4 <- hm_panel(fig_upset(D$setlist, min_size = 1), w = 7.7, h = 4.3,
                title = "Significant-gene overlap (UpSet)")
pC4 <- hm_panel(fig_venn(D$setlist[1:3]), w = 5.1, h = 4.2,
                title = "Significant-gene overlap (Venn)")
pD4 <- hm_panel(fig_lfc_heatmap(D$dfs, n_genes = 40, show_title = FALSE), w = 5.1, h = 4.2,
                title = "Log2 fold-change across contrasts")
p3 <- wrap_plots(
  gg(fig_volcano_grid(D$dfs, n_label = 2, ncol = 2, mode = "publication"),
     "Pairwise differential expression"),
  pB4, pC4, pD4,
  gg(fig_contrast_alluvial(D$dfs, mode = "publication"),
     "Direction of change across contrasts"),
  design = "AAABBB\nCCDDEE", heights = c(1, 1)) + plot_annotation(tag_levels = "A") & tag_theme
save_plate(p3, "figure4", 16, 9.5)
