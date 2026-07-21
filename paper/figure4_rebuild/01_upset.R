## Figure 4b — clean ggplot UpSet, replacing the ComplexHeatmap version.
## Fixes the reported defects: intersection-size bars in green, set-size bars in
## grey (no black), the "Intersection size" axis title placed cleanly with no
## collision against the tick numbers, colorblind-safe, matching the paper theme.
## Reads the cached significant-gene sets (D$setlist) and writes a bare panel PNG.
## Run from the repo root:  Rscript paper/figure4_rebuild/01_upset.R
suppressWarnings(suppressMessages({ library(ggplot2); library(patchwork); library(grid) }))
OUT <- "paper/figure4_rebuild/panels"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
FONT <- "Liberation Sans"
source("R/fig_theme.R")
GREEN <- "#1D9E75"; GREY <- "#9AA3A8"; INK <- "#141414"

D <- readRDS("paper/.panel_cache.rds")
sets <- D$setlist
set_names <- names(sets)
set_sizes <- lengths(sets)
row_levels <- set_names[order(set_sizes, decreasing = TRUE)]   # top row = largest set

universe <- unique(unlist(sets, use.names = FALSE))
M <- vapply(sets, function(s) universe %in% s, logical(length(universe)))
patt <- apply(M, 1L, function(r) paste0(as.integer(r), collapse = ""))
tab  <- sort(table(patt), decreasing = TRUE)
tab  <- tab[names(tab) != paste(rep(0, length(sets)), collapse = "")]
topN <- head(tab, 15L)
combs <- names(topN); sizes <- as.integer(topN); N <- length(combs)

# intersection bar chart (top)
bar_df <- data.frame(x = seq_len(N), size = sizes)
p_bar <- ggplot(bar_df, aes(x = x, y = size)) +
  geom_col(width = 0.62, fill = GREEN) +
  geom_text(aes(label = size), vjust = -0.4, size = 2.3, colour = INK, family = FONT) +
  scale_x_continuous(limits = c(0.5, N + 0.5), expand = c(0, 0)) +
  scale_y_continuous(name = "Intersection size", expand = expansion(mult = c(0, 0.14))) +
  theme_publication() +
  theme(text = element_text(family = FONT),
        axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.y = element_text(size = 12), axis.text.y = element_text(size = 9),
        panel.grid.major.x = element_blank(), plot.margin = margin(4, 6, 2, 4))

# dot matrix (bottom-right)
grid_df <- expand.grid(x = seq_len(N), set = row_levels, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
grid_df$member <- mapply(function(xi, st)
  substr(combs[xi], match(st, set_names), match(st, set_names)) == "1", grid_df$x, grid_df$set)
grid_df$set <- factor(grid_df$set, levels = rev(row_levels))
seg <- do.call(rbind, lapply(seq_len(N), function(xi) {
  members <- row_levels[vapply(row_levels, function(st)
    substr(combs[xi], match(st, set_names), match(st, set_names)) == "1", logical(1))]
  if (length(members) < 2) return(NULL)
  yr <- as.integer(factor(members, levels = rev(row_levels)))
  data.frame(x = xi, y0 = min(yr), y1 = max(yr))
}))
p_mat <- ggplot(grid_df, aes(x = x, y = set)) +
  geom_point(aes(colour = member), size = 3.1) +
  { if (!is.null(seg)) geom_segment(data = seg, aes(x = x, xend = x, y = y0, yend = y1),
                 inherit.aes = FALSE, colour = GREEN, linewidth = 1.0) } +
  scale_colour_manual(values = c(`TRUE` = GREEN, `FALSE` = "#DDE1E3"), guide = "none") +
  scale_x_continuous(limits = c(0.5, N + 0.5), expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = 0.5)) +
  theme_publication() +
  theme(text = element_text(family = FONT), axis.title = element_blank(),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 10, colour = INK),
        panel.grid = element_blank(), plot.margin = margin(2, 6, 4, 4))

# set-size bars (bottom-left, grey, growing leftward)
ss_df <- data.frame(set = factor(row_levels, levels = rev(row_levels)), size = set_sizes[row_levels])
p_ss <- ggplot(ss_df, aes(x = size, y = set)) +
  geom_col(width = 0.62, fill = GREY) +
  scale_x_reverse(name = "Set size", expand = expansion(mult = c(0.06, 0))) +
  scale_y_discrete(expand = expansion(add = 0.5)) +
  theme_publication() +
  theme(text = element_text(family = FONT),
        axis.title.x = element_text(size = 11), axis.text.x = element_text(size = 8),
        axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid.major.y = element_blank(), plot.margin = margin(2, 2, 4, 4))

upset <- p_ss + p_bar + p_mat +
  plot_layout(design = "#B\nAC", widths = c(1.25, 4.2), heights = c(2.1, 1.35))

png(file.path(OUT, "b_upset.png"), width = 6.8, height = 4.6, units = "in", res = 300, bg = "white")
print(upset)
invisible(dev.off()); cat("wrote b_upset (clean ggplot UpSet)\n")
