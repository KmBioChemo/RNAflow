## Figure S1 panel — PROGENy pathway activity (airway, Dex vs control).
## Re-render in the house theme (theme_bare, b=12) so the typography matches
## Figures 2-5, from the cached activity scores. Writes a bare panel PNG.
## Run from the repo root:  Rscript paper/figureS1_rebuild/01_panel.R
suppressWarnings(suppressMessages({ library(ggplot2); library(grid) }))
OUT <- "paper/figureS1_rebuild/panels"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
FONT <- "Liberation Sans"
source("R/fig_theme.R"); source("R/fig_decoupler.R")

theme_bare <- function(b = 12) theme_publication() + theme(
  text = element_text(family = FONT), plot.title = element_blank(),
  axis.title = element_text(size = b), axis.text = element_text(size = b - 2.5),
  legend.title = element_text(size = b - 2, face = "bold"), legend.text = element_text(size = b - 2.5),
  legend.key.size = unit(11, "pt"), plot.margin = margin(3, 6, 3, 4))

D <- readRDS("paper/.activity_cache.rds")
p <- fig_activity_bar(D$act, n = 14, mode = "publication") +
  theme_bare() + theme(legend.position = "right")

# Rendered larger than its final placement (~0.75x) so the apparent font size
# matches the main figures, which place their panels at the same scale.
png(file.path(OUT, "a_activity.png"), width = 9.8, height = 6.3, units = "in", res = 300, bg = "white")
print(p)
invisible(dev.off()); cat("wrote figureS1 a_activity (house theme)\n")
