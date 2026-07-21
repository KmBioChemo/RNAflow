## Figure 2 rework — same treatment as Fig 3.
## egg::ggarrange forces identical plot-area sizes within each row of ggplots
## (so the axes line up); the pheatmap (d) is rendered separately, exactly as
## the GSVA banner is in Fig 3. Panel letters are NOT baked here — they are
## placed uniformly in 02_compose_figure2.py so their size is identical across
## strips regardless of each strip's final scale.
##
## Run from the repo root:  Rscript paper/figure2_rebuild/01_build_panels.R
suppressWarnings(suppressMessages({ library(ggplot2); library(egg); library(grid); library(pheatmap) }))
OUT <- "paper/figure2_rebuild/panels"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
FONT <- "Liberation Sans"
`%||%` <- function(a, b) if (is.null(a)) b else a
source("R/fig_theme.R"); source("R/fig_palettes.R"); source("R/utils_validate.R"); source("R/utils_colors.R")
source("R/fig_volcano.R"); source("R/fig_qc.R"); source("R/fig_enrich.R"); source("R/fig_heatmap.R")

theme_bare <- function(b = 12) theme_publication() + theme(
  text = element_text(family = FONT), plot.title = element_blank(),
  axis.title = element_text(size = b), axis.text = element_text(size = b - 2.5),
  legend.title = element_text(size = b - 2, face = "bold"), legend.text = element_text(size = b - 3),
  legend.position = "bottom", legend.key.size = unit(10, "pt"),
  legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt"),
  plot.margin = margin(3, 6, 3, 6))

D <- readRDS("paper/.panel_cache.rds")

## ---- Row 1 : a volcano | b MA | c p-value histogram --------------------------
a <- fig_volcano(D$a_res, n_label = 8, mode = "publication") + theme_bare() +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2.4)))
b <- fig_ma(D$a_res, mode = "publication") + theme_bare() +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 2.4)))
c <- fig_pval_hist(D$a_res, mode = "publication") + theme_bare()

png(file.path(OUT, "row1.png"), width = 15, height = 4.9, units = "in", res = 300, bg = "white")
egg::ggarrange(a, b, c, nrow = 1, widths = c(1, 1, 1))
invisible(dev.off()); cat("wrote row1 (a,b,c aligned)\n")

## ---- Row 2 : e GSEA dotplot | f ORA bar --------------------------------------
e <- fig_enrich_dot(D$a_gsea, n = 10, mode = "publication") + theme_bare() +
  theme(legend.box = "vertical", legend.box.just = "left",
        legend.spacing.y = unit(1, "pt")) +
  guides(size = guide_legend(title = "Set size", order = 1),
         colour = guide_colourbar(title = "-log10 FDR", order = 2,
           barwidth = unit(2.4, "cm"), barheight = unit(0.28, "cm"),
           title.position = "top", title.hjust = 0))
f <- fig_enrich_bar(D$a_ora, n = 10, mode = "publication") + theme_bare() +
  theme(legend.position = "none")

png(file.path(OUT, "row2.png"), width = 11, height = 4.9, units = "in", res = 300, bg = "white")
egg::ggarrange(e, f, nrow = 1, widths = c(1, 1))
invisible(dev.off()); cat("wrote row2 (e,f aligned)\n")

## ---- d : top differential genes heatmap (pheatmap, rendered separately) ------
hm <- fig_heatmap(D$a_vst, D$a_res, D$am, n_genes = 32, show_colnames = FALSE,
                  direction_annotation = TRUE, show_annotation_names = FALSE,
                  show_title = FALSE)
png(file.path(OUT, "d_heatmap.png"), width = 5.6, height = 4.9, units = "in", res = 300, bg = "white")
grid::grid.draw(hm$gtable)
invisible(dev.off()); cat("wrote d_heatmap\n")
