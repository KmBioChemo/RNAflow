#' Volcano plot functions
#'
#' Static (ggplot2) and interactive (plotly) volcano plots. Both share a
#' common data preparation step and respect the exploration/publication
#' theme mode.
#'
#' @name fig_volcano
NULL

#' Prepare volcano plot data
#'
#' Adds regulation factor, capped -log10(padj), and HTML tooltip text.
#'
#' @param res DE results data.frame (validated)
#' @param lfc_thr log2FoldChange threshold (absolute value)
#' @param padj_thr adjusted p-value threshold
#' @return a tidy data.frame with extra columns: lfc, padj2, nlog, reg, tip
#' @keywords internal
prep_volcano_data <- function(res, lfc_thr, padj_thr) {
  res$lfc   <- as.numeric(res$log2FoldChange)
  res$padj2 <- as.numeric(res$padj)
  res$nlog  <- pmin(-log10(res$padj2 + 1e-300), 320)
  ok <- !is.na(res$padj2) & !is.na(res$lfc)
  reg_chr <- rep("NS", nrow(res))
  reg_chr[ok & res$padj2 < padj_thr & res$lfc >  lfc_thr] <- "Up"
  reg_chr[ok & res$padj2 < padj_thr & res$lfc < -lfc_thr] <- "Down"
  res$reg <- factor(reg_chr, levels = c("Up", "Down", "NS"))
  res$tip <- paste0(
    "<b>", res$gene, "</b><br>",
    "log2FC: ", round(res$lfc, 3), "<br>",
    "FDR: ", formatC(res$padj2, format = "e", digits = 2), "<br>",
    ifelse(res$reg == "Up", "Up-regulated",
           ifelse(res$reg == "Down", "Down-regulated", "Not significant"))
  )
  res
}

#' Static volcano plot (ggplot2)
#'
#' @param res DE results
#' @param lfc_thr log2FoldChange threshold
#' @param padj_thr padj threshold
#' @param n_label number of top genes (by padj) to label with ggrepel
#' @param col_up,col_down,col_ns,col_cut colors for regulation classes and threshold lines
#' @param pt_size point size
#' @param xlab,ylab,title axis and plot labels
#' @param show_title,show_subtitle whether to draw title / subtitle
#' @param leg_pos legend position keyword
#' @param x_min,x_max,y_max optional axis limits
#' @param mode "exploration" or "publication"
#' @return a ggplot2 object
#' @export
fig_volcano <- function(res,
                        lfc_thr = 1, padj_thr = 0.05,
                        n_label = 20,
                        col_up = "#C0392B", col_down = "#2980B9",
                        col_ns = "#BDC3C7", col_cut = "#7F8C8D",
                        pt_size = 1.8,
                        xlab = NULL, ylab = NULL, title = NULL,
                        show_title = TRUE, show_subtitle = TRUE,
                        leg_pos = "Right",
                        x_min = NULL, x_max = NULL, y_max = NULL,
                        mode = c("exploration", "publication")) {

  mode <- match.arg(mode)
  validate_de_results(res)
  df <- prep_volcano_data(res, lfc_thr, padj_thr)

  nu <- sum(df$reg == "Up")
  nd <- sum(df$reg == "Down")
  nn <- sum(df$reg == "NS")
  al <- ceiling(max(abs(df$lfc), na.rm = TRUE) * 1.1)
  xr <- c(if (is_num_scalar(x_min)) x_min else -al,
          if (is_num_scalar(x_max)) x_max else al)

  xlab <- xlab %||% "log2 Fold Change"
  ylab <- ylab %||% "-log10 (adjusted p-value)"

  mt <- if (isTRUE(show_title)) (title %||% "Differential Gene Expression") else NULL
  st <- if (isTRUE(show_subtitle)) sprintf("%d up | %d down | %d NS", nu, nd, nn) else NULL

  # Label top genes
  top <- df[df$reg != "NS", ]
  top <- top[order(top$padj2), ]
  top_genes <- head(top$gene, max(0L, as.integer(n_label)))
  df$lbl <- ifelse(df$gene %in% top_genes, df$gene, NA_character_)

  leg <- switch(
    leg_pos,
    "Top-left"     = ggplot2::theme(legend.position = c(0.01, 0.97),
                                    legend.justification = c(0, 1)),
    "Top-right"    = ggplot2::theme(legend.position = c(0.99, 0.97),
                                    legend.justification = c(1, 1)),
    "Bottom-left"  = ggplot2::theme(legend.position = c(0.01, 0.03),
                                    legend.justification = c(0, 0)),
    "Bottom-right" = ggplot2::theme(legend.position = c(0.99, 0.03),
                                    legend.justification = c(1, 0)),
    "Right"        = ggplot2::theme(legend.position = "right"),
    "None"         = ggplot2::theme(legend.position = "none"),
    ggplot2::theme(legend.position = "right")
  )

  y_scale <- if (is_pos_scalar(y_max)) {
    ggplot2::scale_y_continuous(name = ylab, limits = c(0, y_max),
                                expand = ggplot2::expansion(mult = c(0, 0.03)))
  } else {
    ggplot2::scale_y_continuous(name = ylab,
                                expand = ggplot2::expansion(mult = c(0.01, 0.05)))
  }

  label_size  <- if (mode == "publication") 2.2 else 2.85
  point_scale <- if (mode == "publication") 0.85 else 1

  ggplot2::ggplot(df, ggplot2::aes(x = .data$lfc, y = .data$nlog, color = .data$reg)) +
    ggplot2::geom_point(data = df[df$reg == "NS", ],
                        size = pt_size * 0.70 * point_scale, alpha = 0.35) +
    ggplot2::geom_point(data = df[df$reg != "NS", ],
                        size = pt_size * point_scale, alpha = 0.87) +
    ggplot2::geom_hline(yintercept = -log10(padj_thr),
                        linetype = "dashed", color = col_cut, linewidth = 0.5) +
    ggplot2::geom_vline(xintercept = c(-lfc_thr, lfc_thr),
                        linetype = "dashed", color = col_cut, linewidth = 0.5) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = .data$lbl),
      size = label_size, fontface = "italic", segment.size = 0.28,
      segment.color = "#555", box.padding = 0.38, point.padding = 0.28,
      max.overlaps = 40, na.rm = TRUE
    ) +
    ggplot2::scale_color_manual(
      values = c("Up" = col_up, "Down" = col_down, "NS" = col_ns),
      labels = c("Up"   = sprintf("Up (%d)", nu),
                 "Down" = sprintf("Down (%d)", nd),
                 "NS"   = sprintf("NS (%d)", nn))
    ) +
    ggplot2::scale_x_continuous(
      name = xlab, limits = xr,
      breaks = seq(floor(xr[1]), ceiling(xr[2]),
                   by = max(1, round(diff(xr) / 8)))
    ) +
    y_scale +
    ggplot2::labs(title = mt, subtitle = st, color = NULL) +
    fig_theme(mode) + leg
}

#' Interactive volcano plot (plotly)
#'
#' Same data as [fig_volcano()] but rendered as a plotly figure with hover
#' tooltips and zoom. Always uses exploration-style sizing.
#'
#' @inheritParams fig_volcano
#' @return a plotly object
#' @export
fig_volcano_interactive <- function(res,
                                    lfc_thr = 1, padj_thr = 0.05,
                                    col_up = "#C0392B", col_down = "#2980B9",
                                    col_ns = "#BDC3C7", col_cut = "#7F8C8D",
                                    pt_size = 1.8,
                                    xlab = NULL, ylab = NULL, title = NULL,
                                    show_title = TRUE, show_subtitle = TRUE,
                                    leg_pos = "Right",
                                    x_min = NULL, x_max = NULL, y_max = NULL) {

  validate_de_results(res)
  df <- prep_volcano_data(res, lfc_thr, padj_thr)
  nu <- sum(df$reg == "Up"); nd <- sum(df$reg == "Down"); nn <- sum(df$reg == "NS")
  al <- ceiling(max(abs(df$lfc), na.rm = TRUE) * 1.1)
  xr <- c(if (is_num_scalar(x_min)) x_min else -al,
          if (is_num_scalar(x_max)) x_max else al)
  ytop <- if (is_pos_scalar(y_max)) y_max else NULL

  xlab <- xlab %||% "log2 Fold Change"
  ylab <- ylab %||% "-log10 (adjusted p-value)"

  t_txt <- ""
  if (isTRUE(show_title)) {
    t <- title %||% "Differential Gene Expression"
    if (!nzchar(trimws(t))) t <- "Differential Gene Expression"
    s <- if (isTRUE(show_subtitle)) sprintf("%d up | %d down | %d NS", nu, nd, nn) else ""
    t_txt <- if (nzchar(s)) paste0("<b>", t, "</b><br><sup style='color:#777'>", s, "</sup>")
             else paste0("<b>", t, "</b>")
  } else if (isTRUE(show_subtitle)) {
    t_txt <- sprintf("<span style='color:#777;font-size:.88em;'>%d up | %d down | %d NS</span>", nu, nd, nn)
  }

  cm <- c("Up" = col_up, "Down" = col_down, "NS" = col_ns)
  fig <- plotly::plot_ly(type = "scatter", mode = "markers")
  for (g in c("NS", "Down", "Up")) {
    s <- df[df$reg == g, ]
    if (nrow(s) == 0) next
    fig <- plotly::add_trace(
      fig, x = s$lfc, y = s$nlog,
      name = switch(g,
                    "Up"   = sprintf("Up (%d)", nu),
                    "Down" = sprintf("Down (%d)", nd),
                    sprintf("NS (%d)", nn)),
      text = s$tip, hoverinfo = "text",
      marker = list(color = cm[g],
                    size = if (g == "NS") pt_size * 3.8 else pt_size * 5.2,
                    line = list(width = 0)),
      opacity = if (g == "NS") 0.34 else 0.87
    )
  }
  lpos <- switch(
    leg_pos,
    "Right"        = list(x = 1.01, y = 0.98, xanchor = "left"),
    "Top-right"    = list(x = 0.98, y = 0.98, xanchor = "right"),
    "Top-left"     = list(x = 0.02, y = 0.98, xanchor = "left"),
    "Bottom-right" = list(x = 0.98, y = 0.02, xanchor = "right"),
    "Bottom-left"  = list(x = 0.02, y = 0.02, xanchor = "left"),
    "None"         = list(x = -20, y = 0),
    list(x = 1.01, y = 0.98, xanchor = "left")
  )
  fig %>%
    plotly::layout(
      title = list(text = t_txt, font = list(size = 13), x = 0,
                   xanchor = "left", pad = list(b = 2)),
      xaxis = list(title = xlab, range = xr, zeroline = FALSE,
                   showgrid = FALSE, ticks = "outside"),
      yaxis = list(title = ylab,
                   range = if (!is.null(ytop)) list(0, ytop) else NULL,
                   showgrid = FALSE, ticks = "outside"),
      shapes = list(
        list(type = "line", x0 = -lfc_thr, x1 = -lfc_thr, y0 = 0, y1 = 1,
             yref = "paper", line = list(color = col_cut, dash = "dash", width = 1.2)),
        list(type = "line", x0 =  lfc_thr, x1 =  lfc_thr, y0 = 0, y1 = 1,
             yref = "paper", line = list(color = col_cut, dash = "dash", width = 1.2)),
        list(type = "line", x0 = xr[1] - 1, x1 = xr[2] + 1,
             y0 = -log10(padj_thr), y1 = -log10(padj_thr),
             line = list(color = col_cut, dash = "dash", width = 1.2))
      ),
      legend = c(lpos, list(bgcolor = "rgba(255,255,255,0.9)",
                            bordercolor = "#ddd", borderwidth = 1,
                            font = list(size = 11))),
      hovermode = "closest", paper_bgcolor = "white", plot_bgcolor = "white",
      margin = list(t = if (nzchar(t_txt)) 55 else 22, r = 18, b = 48, l = 58)
    ) %>%
    plotly::config(
      displaylogo = FALSE,
      modeBarButtonsToRemove = c("lasso2d", "select2d", "autoScale2d"),
      toImageButtonOptions = list(format = "png", filename = "volcano_plot", scale = 3)
    )
}
