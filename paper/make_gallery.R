# Render every static figure type RNAflow produces, each as its own standalone
# file (PNG + PDF) in paper/figures/gallery/, using the package's own fig_*
# functions on the bundled datasets. Assemble panels as you wish from these.
# The interactive plotly views (2D/3D PCA, UMAP, linked volcano, enrichment /
# module networks) cannot be exported to a static image and are omitted.
# Run from the package root:  Rscript paper/make_gallery.R
suppressPackageStartupMessages({
  if (!suppressWarnings(require(RNAflow, quietly = TRUE))) {
    devtools::load_all(".", quiet = TRUE)
  }
  library(ggplot2)
})
set.seed(1)
ext    <- function(f) system.file("extdata", f, package = "RNAflow")
outdir <- "paper/figures/gallery"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Save a RNAflow figure (ggplot / pheatmap / ComplexHeatmap) as PNG + PDF.
save_fig <- function(expr, name, w = 6.5, h = 5) {
  obj <- tryCatch(force(expr), error = function(e) {
    message("  SKIP ", name, ": ", conditionMessage(e)); NULL })
  if (is.null(obj)) return(invisible())
  png <- file.path(outdir, paste0(name, ".png"))
  pdf <- file.path(outdir, paste0(name, ".pdf"))
  if (inherits(obj, "ggplot")) {
    ggsave(png, obj, width = w, height = h, dpi = 300, bg = "white")
    ggsave(pdf, obj, width = w, height = h, bg = "white")
  } else if (inherits(obj, "pheatmap")) {
    grDevices::png(png, width = w * 300, height = h * 300, res = 300)
    grid::grid.newpage(); grid::grid.draw(obj$gtable); grDevices::dev.off()
    grDevices::pdf(pdf, width = w, height = h)
    grid::grid.newpage(); grid::grid.draw(obj$gtable); grDevices::dev.off()
  } else if (inherits(obj, c("Heatmap", "HeatmapList"))) {
    grDevices::png(png, width = w * 300, height = h * 300, res = 300)
    ComplexHeatmap::draw(obj); grDevices::dev.off()
    grDevices::pdf(pdf, width = w, height = h)
    ComplexHeatmap::draw(obj); grDevices::dev.off()
  } else if (inherits(obj, c("eulergram", "gTree", "grob", "gList"))) {
    grDevices::png(png, width = w * 300, height = h * 300, res = 300)
    grid::grid.newpage(); grid::grid.draw(obj); grDevices::dev.off()
    grDevices::pdf(pdf, width = w, height = h)
    grid::grid.newpage(); grid::grid.draw(obj); grDevices::dev.off()
  } else { message("  SKIP ", name, ": unhandled class"); return(invisible()) }
  cat("wrote", name, "\n")
}

## ---- inputs ------------------------------------------------------------
a_counts <- read_counts(ext("demo_airway_counts.csv"))
a_meta   <- read_metadata(ext("demo_airway_metadata.csv"), counts_samples = colnames(a_counts))
a_res <- run_deseq2(a_counts, a_meta, design = ~ cell + condition,
                    contrast = c("condition", "Dex", "Control"), shrink = TRUE)
a_sig <- contrast_sig_genes(a_res, padj_thr = 0.05, lfc_thr = 1)
sets  <- get_gene_sets("human", collection = "H")
a_gsea <- run_gsea(a_res, sets, rank_by = "stat")
a_ora  <- run_ora(a_sig, "human", db = "GO", ont = "BP", universe = a_res$gene)

t_counts <- read_counts(ext("demo_tcga_counts.csv"))
t_meta   <- read_metadata(ext("demo_tcga_metadata.csv"), counts_samples = colnames(t_counts))
vst   <- normalize_counts(t_counts, t_meta, "vst")
t_res <- run_deseq2(t_counts, t_meta, design = ~ cancer_type,
                    contrast = c("cancer_type", "LGG", "LUAD"), shrink = TRUE)
gv    <- run_gsva(vst, sets, method = "gsva")

# a small multi-contrast set: a *named list of DE data.frames* (the shape the
# comparison figures expect)
dfs <- list()
for (pr in list(c("LGG","LUAD"), c("KIRC","COAD"), c("BRCA","LUAD"), c("PRAD","THCA"))) {
  dfs[[paste(pr, collapse = " vs ")]] <- run_deseq2(
    t_counts, t_meta, design = ~ cancer_type,
    contrast = c("cancer_type", pr[1], pr[2]), shrink = FALSE)
}
setlist <- contrast_sig_sets(dfs, padj_thr = 0.05, lfc_thr = 1, direction = "either")

datExpr <- wgcna_datexpr(vst, n_genes = 3500)
sft <- wgcna_pick_power(datExpr)
wg  <- run_wgcna(datExpr, power = sft$suggested)
grp <- t_meta$cancer_type[match(rownames(datExpr), t_meta$sample)]

pca <- compute_pca(vst, 500); sc <- merge(pca$scores, t_meta, by = "sample")
OI  <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#7F7F7F")

## ---- differential expression & QC --------------------------------------
save_fig(fig_volcano(a_res, n_label = 12, mode = "publication"), "01_volcano_airway")
save_fig(fig_volcano(t_res, n_label = 12, mode = "publication"), "02_volcano_tcga")
save_fig(fig_pval_hist(a_res, mode = "publication"), "03_pvalue_histogram")
save_fig(fig_ma(a_res, mode = "publication"), "04_ma_plot")
save_fig(fig_lib_sizes(t_counts, t_meta, mode = "publication"), "05_library_sizes", w = 8)
save_fig(fig_sample_cor(vst, t_meta), "06_sample_correlation", w = 7, h = 6)
save_fig(fig_heatmap(vst, t_res, t_meta, n_genes = 30, show_colnames = FALSE),
         "07_de_heatmap", w = 7, h = 6)
save_fig(fig_gene_expression(vst, t_meta, gene = "GFAP", group_by = "cancer_type",
                             style = "raincloud", mode = "publication"),
         "08_gene_expression_raincloud", w = 8)

## ---- sample structure (static rendering of the tool's PCA engine) ------
save_fig(ggplot(sc, aes(PC1, PC2, colour = cancer_type)) +
           geom_point(size = 2.4, alpha = .9) +
           scale_colour_manual(values = OI, name = "Cancer type") +
           labs(x = sprintf("PC1 (%.1f%%)", pca$pct[1]),
                y = sprintf("PC2 (%.1f%%)", pca$pct[2])) + theme_publication(),
         "09_pca_tcga", w = 7)

## ---- functional enrichment ---------------------------------------------
save_fig(fig_enrich_dot(a_gsea, n = 12), "10_gsea_dotplot", w = 8)
save_fig(fig_enrich_bar(a_ora, n = 15, mode = "publication"), "11_ora_barplot", w = 8)
save_fig(fig_gsea_curve(a_res, sets[["HALLMARK_TNFA_SIGNALING_VIA_NFKB"]],
                        title = "HALLMARK TNFA signaling via NFKB", mode = "publication"),
         "12_gsea_running_curve", w = 7)

## ---- co-expression (WGCNA) ---------------------------------------------
save_fig(fig_soft_threshold(sft, mode = "publication"), "13_wgcna_soft_threshold", w = 8)
save_fig(fig_module_sizes(wg, mode = "publication"), "14_wgcna_module_sizes")
save_fig(fig_module_trait(module_trait_cor(wg$MEs, build_traits(t_meta, rownames(datExpr))),
                          mode = "publication"), "15_wgcna_module_trait", w = 8, h = 6)
save_fig(fig_eigengene(wg, "turquoise", groups = grp, mode = "publication"),
         "16_wgcna_eigengene", w = 8)

## ---- per-sample signatures & activity ----------------------------------
save_fig(fig_gsva_heatmap(gv, t_meta, group_by = "cancer_type", n_top = 40),
         "17_gsva_signature_heatmap", w = 9, h = 8)
save_fig(fig_activity_bar(tryCatch(
           run_activity(t_res, get_pathway_network("human"), by = "stat"),
           error = function(e) stop("decoupleR / OmniPath unavailable offline")),
           n = 14, mode = "publication"), "18_pathway_activity")

## ---- multi-contrast comparison -----------------------------------------
save_fig(fig_venn(setlist[1:3]), "19_venn", w = 7)
save_fig(fig_upset(setlist, min_size = 1), "20_upset", w = 9, h = 6)
save_fig(fig_volcano_grid(dfs, n_label = 4, mode = "publication"), "21_volcano_grid", w = 9, h = 7)
save_fig(fig_lfc_heatmap(dfs, n_genes = 40), "22_lfc_heatmap", w = 7, h = 7)
save_fig(fig_contrast_alluvial(dfs, mode = "publication"), "23_direction_alluvial", w = 8)

cat("\nGallery written to ", outdir, "\n")
