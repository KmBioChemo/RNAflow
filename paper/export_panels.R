# Durable panel exporter for the Python montage system (paper/plate/).
# Writes each RNAflow figure panel (Fig 2-4) as a BARE, high-DPI PNG (no title,
# no letter, tight-cropped) into paper/panels/<figure>/. Python (compose.py)
# owns all lettering/titles/layout so the plates are homogeneous. Panels are the
# tool's own outputs (fig_volcano, fig_heatmap, ...), never reconstructions.
# Self-contained: builds paper/.panel_cache.rds on first run, then exports.
# Run:  Rscript paper/export_panels.R
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
ext   <- function(f) system.file("extdata", f, package = "RNAflow")
cache <- "paper/.panel_cache.rds"
# Hallmark sets from the bundled GMT (avoids a msigdbr/msigdbdf download).
read_gmt <- function(f) { l <- strsplit(readLines(f), "\t")
  setNames(lapply(l, function(x) x[-(1:2)]), vapply(l, `[`, "", 1)) }
SETS <- read_gmt("paper/genesets/h.all.v2023.2.Hs.symbols.gmt")
get_gene_sets <- function(...) SETS

## ---- compute once, cache (heavy: DESeq2 / GSVA / WGCNA on the demo data) ----
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

# bare theme: no plot title, tight, uniform text; legend as a compact bottom strip.
# Fonts are rendered large because panels are scaled down when composed into the
# plate -- larger source text stays legible in the final montage.
theme_bare <- function(b = 14) theme_publication() + theme(
  text = element_text(family = FONT),
  plot.title = element_blank(),
  axis.title.x = element_text(size = b), axis.title.y = element_text(size = b),
  axis.text = element_text(size = b - 2.5),
  legend.title = element_text(size = b - 1.5, face = "bold"),
  legend.text = element_text(size = b - 2.5),
  legend.position = "bottom", legend.key.size = unit(12, "pt"),
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
# top row (a,b,c): square-ish. bottom row (d,e,f) share one aspect so their
# heights match the heatmap when composed.
save_gg(fig_volcano(D$a_res, n_label = 8, mode = "publication") + theme_bare(), dir2, "a_volcano", 5.6, 5.0)
save_gg(fig_ma(D$a_res, mode = "publication") + theme_bare(),                    dir2, "b_ma", 5.6, 5.0)
save_gg(fig_pval_hist(D$a_res, mode = "publication") + theme_bare(),             dir2, "c_pval", 5.6, 5.0)
save_hm(fig_heatmap(D$a_vst, D$a_res, D$am, n_genes = 32, show_colnames = FALSE,
                    direction_annotation = TRUE, show_annotation_names = FALSE,
                    show_title = FALSE),                                          dir2, "d_heatmap", 5.6, 4.9)
# GSEA as a bubble/dot plot (NES on x, set size = dot size, colour = -log10 FDR);
# legends stacked vertically on the right
save_gg(fig_enrich_dot(D$a_gsea, n = 10, mode = "publication") + theme_bare() +
          theme(legend.position = "right", legend.box = "vertical",
                legend.direction = "vertical"),
        dir2, "e_gsea", 5.6, 4.9)
save_gg(fig_enrich_bar(D$a_ora, n = 10, mode = "publication") + theme_bare() +
          theme(legend.position = "none"),                                       dir2, "f_ora", 5.6, 4.9)

## ===== figure 3 : molecular landscape & co-expression (TCGA) ================
dir3 <- "paper/panels/figure3"; dir.create(dir3, showWarnings = FALSE, recursive = TRUE)
cat("figure3:\n")
# GSVA hero: rendered TALL (aspect ~ the tall hero cell) so the heatmap fills the
# full height of its panel instead of letterboxing.
save_hm(fig_gsva_heatmap(D$gv, D$tm, group_by = "cancer_type", n_top = 45, title = "",
                         show_annotation_names = FALSE),
        dir3, "a_gsva", 6.4, 9.3)
sc <- merge(D$pca$scores, D$tm, by = "sample")
p_pca <- ggplot(sc, aes(PC1, PC2, colour = cancer_type)) +
  geom_point(size = 2.4, alpha = .92, stroke = 0) +
  scale_colour_manual(values = OI, name = "Cancer type") +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 2.6))) +
  labs(x = sprintf("PC1 (%.1f%%)", D$pca$pct[1]), y = sprintf("PC2 (%.1f%%)", D$pca$pct[2])) +
  theme_bare()
save_gg(p_pca, dir3, "b_pca", 4.6, 4.7)
# soft-threshold: two facets -> render wide enough that the power labels breathe
save_gg(fig_soft_threshold(D$sft, mode = "publication") + theme_bare(), dir3, "c_soft", 5.0, 4.4)
save_gg(fig_module_trait(D$mt, mode = "publication") + theme_bare() +
          theme(axis.text.x = element_text(angle = 40, hjust = 1, size = 9),
                legend.position = "right", legend.direction = "vertical"),
        dir3, "d_modtrait", 4.8, 4.6)
# module enrichment: dense panel -> smaller base font, rotated coloured module
# labels, compact legends, rendered wide so nothing is clipped
p2E <- if (!is.null(D$mod_enrich))
  fig_module_enrichment(D$mod_enrich, max_terms = 14, mode = "publication") + theme_bare(11) +
    theme(axis.text.x = ggtext::element_markdown(size = 9, angle = 30, hjust = 1),
          axis.text.y = element_text(size = 9),
          legend.position = "right", legend.box = "vertical", legend.direction = "vertical",
          legend.text = element_text(size = 8), legend.title = element_text(size = 9, face = "bold"),
          legend.key.size = unit(9, "pt")) else
  fig_module_sizes(D$wg, mode = "publication") + theme_bare()
save_gg(p2E, dir3, "e_modenrich", 6.0, 4.8)

## ===== figure 4 : multi-contrast comparison (TCGA) ==========================
dir4 <- "paper/panels/figure4"; dir.create(dir4, showWarnings = FALSE, recursive = TRUE)
cat("figure4:\n")
save_gg(fig_volcano_grid(D$dfs, n_label = 2, ncol = 2, mode = "publication") + theme_bare(),
        dir4, "a_volcgrid", 7.6, 4.6)
# UpSet: wider set-size bars so they breathe next to the intersection matrix
save_hm(fig_upset(D$setlist, min_size = 1, set_size_width = 3.2), dir4, "b_upset", 6.8, 4.6)
# Venn: render on a larger canvas -> more air around the circles
save_hm(fig_venn(D$setlist[1:3]), dir4, "c_venn", 5.2, 5.0)
save_hm(fig_lfc_heatmap(D$dfs, n_genes = 40, show_title = FALSE), dir4, "d_lfc", 4.9, 4.8)
save_gg(fig_contrast_alluvial(D$dfs, mode = "publication") + theme_bare() +
          theme(axis.text.x = element_text(angle = 20, hjust = 1, size = 10)),
        dir4, "e_alluvial", 5.4, 4.8)

cat("done: panels in paper/panels/\n")
