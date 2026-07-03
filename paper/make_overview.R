# Graphical-abstract overview (Figure 1): a feature-map schematic of the whole
# RNAflow workflow (input -> all analysis modules -> reproducible outputs) plus a
# strip of real output thumbnails, so a reader sees every capability at a glance.
# Run from the package root:  Rscript paper/make_overview.R
suppressPackageStartupMessages({ library(ggplot2); library(cowplot) })

## ---- schematic --------------------------------------------------------
w <- 17; h <- 11
chips <- data.frame(
  label = c("Differential\nexpression", "Quality\ncontrol", "Sample overview\n(PCA / UMAP / 3D)",
            "Functional\nenrichment", "Co-expression\n(WGCNA)", "Activity\n(decoupleR)",
            "Per-sample\nsignatures (GSVA)", "Multi-contrast\ncomparison", "AI-assisted\ninterpretation"),
  col  = rep(c("#e9f1fb", "#e7f4ee", "#fdefe6"), each = 3),
  bord = rep(c("#4a90d9", "#1D9E75", "#D55E00"), each = 3),
  cx = rep(c(31, 50, 69), times = 3),
  cy = rep(c(69, 53.5, 38), each = 3))

sch <- ggplot() +
  # input
  annotate("rect", xmin = 3, xmax = 17, ymin = 45, ymax = 62, fill = "#eceff3",
           colour = "#8a94a3", linewidth = .6) +
  annotate("text", x = 10, y = 53.5, label = "Counts +\nmetadata", size = 3.5,
           fontface = "bold", lineheight = .9) +
  annotate("text", x = 10, y = 41, label = "human · mouse · rat", size = 2.5, colour = "#7a828e") +
  # analyses container + title
  annotate("rect", xmin = 21.5, xmax = 78.5, ymin = 29, ymax = 79, fill = NA,
           colour = "#c7ced8", linewidth = .5) +
  annotate("text", x = 50, y = 82.5, label = "ANALYSES  (one loaded dataset)", size = 3.6,
           fontface = "bold", colour = "#5a6472") +
  geom_rect(data = chips, aes(xmin = cx - w/2, xmax = cx + w/2, ymin = cy - h/2, ymax = cy + h/2,
                              fill = I(col), colour = I(bord)), linewidth = .6) +
  geom_text(data = chips, aes(cx, cy, label = label), size = 2.95, lineheight = .85,
            fontface = "bold", colour = "#22303a") +
  # outputs
  annotate("rect", xmin = 83, xmax = 99, ymin = 36, ymax = 71, fill = "#eceff3",
           colour = "#8a94a3", linewidth = .6) +
  annotate("text", x = 91, y = 66, label = "Reproducible\nexport", size = 3.3,
           fontface = "bold", lineheight = .9) +
  annotate("text", x = 91, y = 51,
           label = "R script\nMethods paragraph\nself-contained report\nsession (.rds)\nDocker image",
           size = 2.75, lineheight = 1.15) +
  # arrows
  annotate("segment", x = 17.3, xend = 21.2, y = 53.5, yend = 53.5, colour = "#5a6472",
           linewidth = .8, arrow = arrow(length = unit(2.4, "mm"), type = "closed")) +
  annotate("segment", x = 78.8, xend = 82.7, y = 53.5, yend = 53.5, colour = "#5a6472",
           linewidth = .8, arrow = arrow(length = unit(2.4, "mm"), type = "closed")) +
  coord_cartesian(xlim = c(0, 100), ylim = c(25, 87), expand = FALSE) +
  theme_void()

## ---- thumbnail strip (real outputs) -----------------------------------
gal <- "paper/figures/gallery"
thumbs <- c("01_volcano_airway", "09_pca_tcga", "07_de_heatmap", "10_gsea_dotplot",
            "15_wgcna_module_trait", "17_gsva_signature_heatmap", "20_upset", "23_direction_alluvial")
labs <- c("Volcano", "PCA", "Heatmap", "Enrichment",
          "WGCNA module–trait", "GSVA signatures", "UpSet", "Direction alluvial")
tp <- lapply(seq_along(thumbs), function(i) {
  p <- file.path(gal, paste0(thumbs[i], ".png"))
  ggdraw() +
    draw_image(p, x = 0, y = 0.11, width = 1, height = 0.89) +
    draw_label(labs[i], y = 0.045, size = 8, fontface = "bold", colour = "#2b2f36")
})
strip <- plot_grid(plotlist = tp, nrow = 1)

## ---- compose ----------------------------------------------------------
fig <- plot_grid(sch, strip, ncol = 1, rel_heights = c(2.15, 1))
ggsave("paper/figures/figure1.png", fig, width = 14, height = 8, dpi = 300, bg = "white")
ggsave("paper/figures/figure1.pdf", fig, width = 14, height = 8, bg = "white")
cat("wrote figure1\n")
