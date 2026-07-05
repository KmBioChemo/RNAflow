# Export each Figure 2 plot as its own file, all at IDENTICAL dimensions, for
# manual montage. Panels are the RNAflow tool's own outputs (no reconstruction).
# Run:  Rscript paper/make_fig2_panels.R   (needs paper/.panel_cache.rds)
suppressPackageStartupMessages({
  if (!suppressWarnings(require(RNAflow, quietly = TRUE))) devtools::load_all(".", quiet = TRUE)
  library(ggplot2); library(patchwork); library(ragg); library(png)
})
set.seed(1)
FONT <- { pref <- c("Helvetica","Arial","Liberation Sans","DejaVu Sans")
  fams <- tryCatch(systemfonts::system_fonts()$family, error=function(e) character(0))
  hit <- pref[pref %in% fams]; if (length(hit)) hit[1] else "sans" }
D <- readRDS("paper/.panel_cache.rds")
if (exists("get_gene_sets")) sets <- get_gene_sets("human", collection = "H") else
  sets <- { l <- strsplit(readLines("paper/genesets/h.all.v2023.2.Hs.symbols.gmt"), "\t")
            setNames(lapply(l, function(x) x[-(1:2)]), vapply(l, `[`, "", 1)) }

## ---- uniform panel geometry -----------------------------------------------
PW <- 5.5; PH <- 5.0                      # every panel is exactly this size (inches)
OUT <- "paper/figures/fig2_panels"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

theme_2026 <- function(b = 11.5) theme_publication() + theme(
  text = element_text(family = FONT),
  plot.title = element_text(size = b - 0.5, face = "bold", margin = margin(b = 2)),
  axis.title.x = element_text(size = b - 2, margin = margin(t = 1)),
  axis.title.y = element_text(size = b - 2, margin = margin(r = 1)),
  axis.text = element_text(size = b - 3.5),
  legend.title = element_text(size = b - 3, face = "bold"),
  legend.text = element_text(size = b - 4),
  legend.key.size = unit(9, "pt"), legend.margin = margin(0, 0, 0, 0),
  legend.box.spacing = unit(3, "pt"), plot.margin = margin(4, 5, 4, 5))
gg   <- function(p, title) p + labs(title = title) + theme_2026()
gg_b <- function(p, title) gg(p, title) +
  theme(legend.position = "bottom", legend.title = element_text(size = 8, face = "bold"),
        legend.text = element_text(size = 7.5), legend.box.spacing = unit(2, "pt"))

.trim_white <- function(img, tol = 0.985, pad = 5) {
  d <- dim(img); if (length(d) == 2) img <- array(img, c(d, 1))
  ch <- min(dim(img)[3], 3)
  ink <- Reduce(`|`, lapply(seq_len(ch), function(k) img[, , k] < tol))
  rows <- which(rowSums(ink) > 0); cols <- which(colSums(ink) > 0)
  if (!length(rows) || !length(cols)) return(img)
  img[max(1,min(rows)-pad):min(dim(img)[1],max(rows)+pad),
      max(1,min(cols)-pad):min(dim(img)[2],max(cols)+pad), , drop = FALSE]
}
hm_panel <- function(obj, w, h, title, dpi = 320) {
  f <- tempfile(fileext = ".png")
  ragg::agg_png(f, width = w, height = h, units = "in", res = dpi, background = "white")
  if (inherits(obj, "pheatmap")) grid::grid.draw(obj$gtable)
  else if (inherits(obj, c("Heatmap","HeatmapList"))) ComplexHeatmap::draw(obj)
  else grid::grid.draw(obj)
  grDevices::dev.off()
  img <- .trim_white(png::readPNG(f)); W <- dim(img)[2]; H <- dim(img)[1]
  ggplot() +
    annotation_custom(grid::rasterGrob(img, width = unit(1,"npc"), height = unit(1,"npc")), 0, W, 0, H) +
    coord_fixed(1, xlim = c(0, W), ylim = c(0, H), expand = FALSE, clip = "off") +
    theme_void() + theme(plot.margin = margin(1, 1, 1, 1)) +
    ggtitle(title) + theme(plot.title = element_text(size = 11, face = "bold",
      family = FONT, hjust = 0, margin = margin(b = 2)))
}

## ---- the six panels (same order/content as figure 2) ----------------------
panels <- list(
  A_differential_expression = gg_b(fig_volcano(D$a_res, n_label = 8, mode = "publication"),
                                    "Differential expression"),
  B_ma_plot                 = gg_b(fig_ma(D$a_res, mode = "publication"), "MA plot"),
  C_pvalue_histogram        = gg(fig_pval_hist(D$a_res, mode = "publication"), "P-value histogram"),
  D_top_genes               = hm_panel(fig_heatmap(D$a_vst, D$a_res, D$am, n_genes = 32,
                                          show_colnames = FALSE, direction_annotation = TRUE,
                                          show_title = FALSE), w = PW, h = PH - 0.4,
                                        title = "Top differential genes"),
  E_gsea                    = gg_b(fig_gsea_ridge(D$a_res, sets, D$a_gsea, n = 10, mode = "publication"),
                                   "Gene-set enrichment (GSEA)"),
  F_ora                     = gg(fig_enrich_bar(D$a_ora, n = 10, mode = "publication"),
                                 "Over-representation (GO BP)") + theme(legend.position = "none"))

for (nm in names(panels)) {
  ggsave(file.path(OUT, paste0(nm, ".png")), panels[[nm]], width = PW, height = PH,
         dpi = 300, bg = "white", device = ragg::agg_png)
  ggsave(file.path(OUT, paste0(nm, ".pdf")), panels[[nm]], width = PW, height = PH,
         bg = "white", device = grDevices::cairo_pdf)
  cat("wrote", nm, sprintf("(%.1f x %.1f in)\n", PW, PH))
}
cat("all panels in", OUT, "\n")
