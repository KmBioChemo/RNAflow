# Durable panel exporter for the Python montage system (paper/plate/).
# Writes each RNAflow figure panel as a BARE, high-DPI PNG (no title, no letter,
# tight-cropped) into paper/panels/<figure>/. Python (compose.py) owns all
# lettering/titles/layout so the plates are homogeneous. Panels are the tool's
# own outputs (fig_volcano, fig_heatmap, ...), never reconstructions.
# Run:  Rscript paper/export_panels.R            (needs paper/.panel_cache.rds)
suppressPackageStartupMessages({
  if (!suppressWarnings(require(RNAflow, quietly = TRUE))) devtools::load_all(".", quiet = TRUE)
  library(ggplot2); library(ragg); library(png)
})
set.seed(1)
FONT <- { pref <- c("Helvetica","Arial","Liberation Sans","DejaVu Sans")
  fams <- tryCatch(systemfonts::system_fonts()$family, error=function(e) character(0))
  hit <- pref[pref %in% fams]; if (length(hit)) hit[1] else "sans" }
DPI <- 400
OI <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#7F7F7F")
D <- readRDS("paper/.panel_cache.rds")
if (exists("get_gene_sets")) sets <- get_gene_sets("human", collection = "H") else
  sets <- { l <- strsplit(readLines("paper/genesets/h.all.v2023.2.Hs.symbols.gmt"), "\t")
            setNames(lapply(l, function(x) x[-(1:2)]), vapply(l, `[`, "", 1)) }

# bare theme: no plot title, tight, uniform text; legend as a compact bottom strip
theme_bare <- function(b = 11) theme_publication() + theme(
  text = element_text(family = FONT),
  plot.title = element_blank(),
  axis.title.x = element_text(size = b - 1.5), axis.title.y = element_text(size = b - 1.5),
  axis.text = element_text(size = b - 3),
  legend.title = element_text(size = b - 2.5, face = "bold"),
  legend.text = element_text(size = b - 3),
  legend.position = "bottom", legend.key.size = unit(10, "pt"),
  legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
  plot.margin = margin(3, 4, 3, 4))

# crop uniform white border so the panel is tight
.trim <- function(img, tol = 0.985, pad = 4) {
  d <- dim(img); if (length(d) == 2) img <- array(img, c(d, 1))
  ch <- min(dim(img)[3], 3)
  ink <- Reduce(`|`, lapply(seq_len(ch), function(k) img[, , k] < tol))
  r <- which(rowSums(ink) > 0); c <- which(colSums(ink) > 0)
  if (!length(r) || !length(c)) return(img)
  img[max(1,min(r)-pad):min(dim(img)[1],max(r)+pad),
      max(1,min(c)-pad):min(dim(img)[2],max(c)+pad), , drop = FALSE]
}
save_gg <- function(p, dir, name, w, h) {
  f <- file.path(dir, paste0(name, ".png"))
  ggsave(f, p, width = w, height = h, dpi = DPI, bg = "white", device = ragg::agg_png)
  img <- .trim(png::readPNG(f)); png::writePNG(img, f)          # tight-crop in place
  cat("  ", name, sprintf("(%dx%d px)\n", dim(img)[2], dim(img)[1]))
}
save_hm <- function(obj, dir, name, w, h) {                     # pheatmap/ComplexHeatmap/grob
  f <- file.path(dir, paste0(name, ".png"))
  ragg::agg_png(f, width = w, height = h, units = "in", res = DPI, background = "white")
  if (inherits(obj, "pheatmap")) grid::grid.draw(obj$gtable)
  else if (inherits(obj, c("Heatmap","HeatmapList"))) ComplexHeatmap::draw(obj)
  else grid::grid.draw(obj)
  grDevices::dev.off()
  img <- .trim(png::readPNG(f)); png::writePNG(img, f)
  cat("  ", name, sprintf("(%dx%d px)\n", dim(img)[2], dim(img)[1]))
}

## ===== figure 2 : differential expression & enrichment (airway) =============
dir2 <- "paper/panels/figure2"; dir.create(dir2, showWarnings = FALSE, recursive = TRUE)
cat("figure2:\n")
save_gg(fig_volcano(D$a_res, n_label = 8, mode = "publication") + theme_bare(), dir2, "a_volcano", 5.6, 5.0)
save_gg(fig_ma(D$a_res, mode = "publication") + theme_bare(),                    dir2, "b_ma", 5.6, 5.0)
save_gg(fig_pval_hist(D$a_res, mode = "publication") + theme_bare(),             dir2, "c_pval", 5.6, 5.0)
save_hm(fig_heatmap(D$a_vst, D$a_res, D$am, n_genes = 32, show_colnames = FALSE,
                    direction_annotation = TRUE, show_annotation_names = FALSE,
                    show_title = FALSE),                                          dir2, "d_heatmap", 5.6, 5.0)
save_gg(fig_gsea_ridge(D$a_res, sets, D$a_gsea, n = 10, mode = "publication") + theme_bare(),
        dir2, "e_gsea", 5.6, 5.0)
save_gg(fig_enrich_bar(D$a_ora, n = 10, mode = "publication") + theme_bare() +
          theme(legend.position = "none"),                                       dir2, "f_ora", 5.6, 5.0)

## ===== figure 3 : molecular landscape & co-expression (TCGA) ================
dir3 <- "paper/panels/figure3"; dir.create(dir3, showWarnings = FALSE, recursive = TRUE)
cat("figure3:\n")
# GSVA hero rendered near-square to fill the 2x2 hero cell (many rows + samples)
save_hm(fig_gsva_heatmap(D$gv, D$tm, group_by = "cancer_type", n_top = 45, title = "",
                         show_annotation_names = FALSE),
        dir3, "a_gsva", 7.8, 7.3)
sc <- merge(D$pca$scores, D$tm, by = "sample")
p_pca <- ggplot(sc, aes(PC1, PC2, colour = cancer_type)) +
  geom_point(size = 2.4, alpha = .92, stroke = 0) +
  scale_colour_manual(values = OI, name = "Cancer type") +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 2.6))) +
  labs(x = sprintf("PC1 (%.1f%%)", D$pca$pct[1]), y = sprintf("PC2 (%.1f%%)", D$pca$pct[2])) +
  theme_bare()
save_gg(p_pca, dir3, "b_pca", 4.7, 4.7)
save_gg(fig_soft_threshold(D$sft, mode = "publication") + theme_bare(), dir3, "c_soft", 4.7, 4.7)
save_gg(fig_module_trait(D$mt, mode = "publication") + theme_bare() +
          theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 6.5)),
        dir3, "d_modtrait", 4.7, 4.7)
p2E <- if (!is.null(D$mod_enrich))
  fig_module_enrichment(D$mod_enrich, max_terms = 15, mode = "publication") + theme_bare() +
    theme(axis.text.x = ggtext::element_markdown(size = 7)) else
  fig_module_sizes(D$wg, mode = "publication") + theme_bare()
save_gg(p2E, dir3, "e_modenrich", 4.7, 4.7)

## ===== figure 4 : multi-contrast comparison (TCGA) ==========================
dir4 <- "paper/panels/figure4"; dir.create(dir4, showWarnings = FALSE, recursive = TRUE)
cat("figure4:\n")
save_gg(fig_volcano_grid(D$dfs, n_label = 2, ncol = 2, mode = "publication") + theme_bare(),
        dir4, "a_volcgrid", 7.6, 4.6)
save_hm(fig_upset(D$setlist, min_size = 1), dir4, "b_upset", 7.6, 4.4)
save_hm(fig_venn(D$setlist[1:3]), dir4, "c_venn", 4.9, 4.5)
save_hm(fig_lfc_heatmap(D$dfs, n_genes = 40, show_title = FALSE), dir4, "d_lfc", 4.9, 4.6)
save_gg(fig_contrast_alluvial(D$dfs, mode = "publication") + theme_bare(), dir4, "e_alluvial", 5.4, 4.6)

cat("done: panels in paper/panels/\n")
