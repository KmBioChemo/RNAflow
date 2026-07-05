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
cat("done: panels in paper/panels/\n")
