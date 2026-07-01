#' WGCNA figures
#'
#' Plots for the co-expression module: soft-threshold diagnostics, the
#' module-trait correlation heatmap, module sizes, and module eigengene
#' profiles. Pure ggplot2 functions respecting the theme mode.
#'
#' @name fig_wgcna
NULL

#' Soft-threshold diagnostic plot
#'
#' Scale-free topology fit (signed R^2) and mean connectivity vs. power, with
#' the R^2 target line and the suggested power highlighted.
#'
#' @param sft the list returned by [wgcna_pick_power()]
#' @param mode "exploration" or "publication"
#' @return a ggplot2 object (two facets)
#' @export
fig_soft_threshold <- function(sft, mode = c("exploration", "publication")) {
  mode <- match.arg(mode)
  fi <- sft$fit_indices
  rsq <- -sign(fi$slope) * fi$SFT.R.sq
  long <- rbind(
    data.frame(power = fi$Power, value = rsq,
               metric = "Scale-free fit (signed R2)"),
    data.frame(power = fi$Power, value = fi$mean.k.,
               metric = "Mean connectivity")
  )
  long$metric <- factor(long$metric,
                        levels = c("Scale-free fit (signed R2)",
                                   "Mean connectivity"))
  hl <- data.frame(metric = factor("Scale-free fit (signed R2)",
                                    levels = levels(long$metric)),
                   yint = sft$rsq_cut)

  ggplot2::ggplot(long, ggplot2::aes(.data$power, .data$value)) +
    ggplot2::geom_hline(data = hl, ggplot2::aes(yintercept = .data$yint),
                        linetype = "dashed", color = "#C0392B", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = sft$suggested, linetype = "dotted",
                        color = "#1D9E75", linewidth = 0.5) +
    ggplot2::geom_line(color = "#95A5A6", linewidth = 0.4) +
    ggplot2::geom_text(ggplot2::aes(label = .data$power), size = 2.8,
                       color = "#2C3E50") +
    ggplot2::facet_wrap(~ .data$metric, scales = "free_y") +
    ggplot2::labs(x = "Soft-threshold power", y = NULL,
                  title = sprintf("Suggested power: %s", sft$suggested)) +
    fig_theme(mode)
}

#' Module-trait correlation heatmap
#'
#' @param mt the list returned by [module_trait_cor()]
#' @param mode "exploration" or "publication"
#' @param text_size size of the in-cell r / p labels
#' @return a ggplot2 object
#' @export
fig_module_trait <- function(mt, mode = c("exploration", "publication"),
                             text_size = 2.6) {
  mode <- match.arg(mode)
  cmat <- mt$cor; pmat <- mt$p
  df <- expand.grid(module = rownames(cmat), trait = colnames(cmat),
                    stringsAsFactors = FALSE)
  df$r <- as.vector(cmat)
  df$p <- as.vector(pmat)
  df$module <- sub("^ME", "", df$module)
  df$label <- sprintf("%.2f\n(%.0e)", df$r, df$p)
  df$txt_col <- ifelse(abs(df$r) > 0.6, "white", "black")
  df$module <- factor(df$module, levels = rev(sub("^ME", "", rownames(cmat))))
  df$trait  <- factor(df$trait, levels = colnames(cmat))

  ggplot2::ggplot(df, ggplot2::aes(.data$trait, .data$module, fill = .data$r)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label, color = .data$txt_col),
                       size = text_size, lineheight = 0.85) +
    ggplot2::scale_color_identity() +
    ggplot2::scale_fill_gradient2(low = "#2980B9", mid = "white",
                                  high = "#C0392B", midpoint = 0,
                                  limits = c(-1, 1), name = "Correlation") +
    ggplot2::labs(x = NULL, y = "Module") +
    fig_theme(mode) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Module size bar chart
#'
#' Bars colored by the module's own WGCNA color.
#'
#' @param wg the list returned by [run_wgcna()]
#' @param include_grey include the unassigned grey module
#' @param mode "exploration" or "publication"
#' @return a ggplot2 object
#' @export
fig_module_sizes <- function(wg, include_grey = FALSE,
                             mode = c("exploration", "publication")) {
  mode <- match.arg(mode)
  ms <- module_summary(wg)
  if (!include_grey) ms <- ms[ms$module != "grey", , drop = FALSE]
  if (nrow(ms) == 0) stop("No modules to plot.", call. = FALSE)
  ms$module <- factor(ms$module, levels = ms$module[order(ms$n_genes)])

  ggplot2::ggplot(ms, ggplot2::aes(.data$n_genes, .data$module,
                                   fill = .data$module)) +
    ggplot2::geom_col(width = 0.7, color = "#444444", linewidth = 0.25) +
    ggplot2::geom_text(ggplot2::aes(label = .data$n_genes),
                       hjust = -0.2, size = 2.8) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(x = "Genes in module", y = NULL) +
    fig_theme(mode)
}

#' Module eigengene profile
#'
#' Per-sample eigengene value for one module, optionally grouped/colored by a
#' sample annotation.
#'
#' @param wg the list returned by [run_wgcna()]
#' @param module module color
#' @param groups optional vector of group labels aligned to the samples
#'   (`rownames(wg$MEs)`); if named, it is reordered to match
#' @param mode "exploration" or "publication"
#' @return a ggplot2 object
#' @export
fig_eigengene <- function(wg, module, groups = NULL,
                          mode = c("exploration", "publication")) {
  mode <- match.arg(mode)
  me_col <- paste0("ME", module)
  if (!me_col %in% colnames(wg$MEs)) {
    stop("Module '", module, "' has no eigengene.", call. = FALSE)
  }
  df <- data.frame(sample = rownames(wg$MEs),
                   eigengene = wg$MEs[[me_col]],
                   stringsAsFactors = FALSE)
  if (!is.null(groups)) {
    g <- if (!is.null(names(groups))) groups[df$sample] else groups
    df$group <- as.character(g)
    df <- df[order(df$group, df$eigengene), ]
  } else {
    df <- df[order(df$eigengene), ]
  }
  df$sample <- factor(df$sample, levels = df$sample)

  p <- ggplot2::ggplot(df, ggplot2::aes(.data$sample, .data$eigengene))
  if (!is.null(groups)) {
    p <- p + ggplot2::geom_col(ggplot2::aes(fill = .data$group), width = 0.8) +
      ggplot2::labs(fill = NULL)
  } else {
    p <- p + ggplot2::geom_col(fill = module, color = "#444444",
                               linewidth = 0.2, width = 0.8)
  }
  p +
    ggplot2::labs(x = NULL, y = "Eigengene",
                  title = sprintf("Module %s eigengene", module)) +
    fig_theme(mode) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1,
                                                       vjust = 0.5, size = 6))
}
