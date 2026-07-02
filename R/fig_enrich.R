#' Functional enrichment figures
#'
#' Plots for GSEA and ORA results: dotplot, lollipop bar, and the GSEA
#' running-enrichment curve. Pure ggplot2 functions that respect the
#' exploration / publication theme mode.
#'
#' @name fig_enrich
#' @keywords internal
NULL

#' Tidy enrichment label
#'
#' Strips collection prefixes (e.g. "HALLMARK_") and turns underscores into
#' spaces for readable axis labels.
#'
#' @param x character vector of term names
#' @keywords internal
clean_term <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  trimws(x)
}

#' Wrap long labels onto multiple lines
#' @keywords internal
wrap_label <- function(x, width = 42) {
  vapply(x, function(s) paste(strwrap(s, width = width), collapse = "\n"),
         character(1), USE.NAMES = FALSE)
}

#' Normalize a GSEA or ORA table to a common plotting frame
#'
#' @param df a data.frame from [run_gsea()] or [run_ora()]
#' @param n keep the top `n` terms by padj
#' @return a data.frame with columns: term, score, count, padj, direction,
#'   is_gsea
#' @keywords internal
enrich_plot_df <- function(df, n = 20) {
  if (!is.data.frame(df) || nrow(df) == 0) {
    stop("No enrichment results to plot.", call. = FALSE)
  }
  is_gsea <- "NES" %in% colnames(df)
  if (is_gsea) {
    out <- data.frame(
      term = as.character(df$pathway),
      score = as.numeric(df$NES),
      count = as.numeric(df$size %||% NA),
      padj = as.numeric(df$padj),
      direction = ifelse(as.numeric(df$NES) >= 0, "Up", "Down"),
      is_gsea = TRUE, stringsAsFactors = FALSE)
  } else {
    ratio <- vapply(strsplit(as.character(df$GeneRatio %||% ""), "/"),
                    function(p) if (length(p) == 2)
                      as.numeric(p[1]) / as.numeric(p[2]) else NA_real_,
                    numeric(1))
    out <- data.frame(
      term = as.character(df$Description),
      score = if (all(is.na(ratio))) as.numeric(df$Count) else ratio,
      count = as.numeric(df$Count),
      padj = as.numeric(df$padj),
      direction = "Enriched",
      is_gsea = FALSE, stringsAsFactors = FALSE)
  }
  out <- out[order(out$padj), , drop = FALSE]
  utils::head(out, max(1L, as.integer(n)))
}

#' Enrichment dotplot
#'
#' @param df a data.frame from [run_gsea()] or [run_ora()]
#' @param n number of top terms (by padj) to show
#' @param col_up,col_down colors for positive / negative NES (GSEA only)
#' @param mode "exploration" or "publication"
#' @return a ggplot2 object
#' @export
fig_enrich_dot <- function(df, n = 20,
                           col_up = "#C0392B", col_down = "#2980B9",
                           mode = c("exploration", "publication")) {
  mode <- match.arg(mode)
  d <- enrich_plot_df(df, n)
  d$lbl <- clean_term(d$term)
  d$lbl <- factor(d$lbl, levels = d$lbl[order(d$score)])
  xlab <- if (d$is_gsea[1]) "Normalized enrichment score (NES)" else "Gene ratio"

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$score, y = .data$lbl)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = .data$score,
                   y = .data$lbl, yend = .data$lbl),
      color = "#CCCCCC", linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(size = .data$count,
                                     color = -log10(.data$padj))) +
    ggplot2::scale_color_viridis_c(option = "C", end = 0.9,
                                   name = expression(-log[10] ~ FDR)) +
    ggplot2::scale_size_continuous(name = "Set size", range = c(2, 7)) +
    ggplot2::scale_y_discrete(labels = function(x) wrap_label(x)) +
    ggplot2::labs(x = xlab, y = NULL) +
    fig_theme(mode)

  if (d$is_gsea[1]) {
    p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                                 color = "#7F8C8D", linewidth = 0.4)
  }
  p
}

#' Enrichment lollipop bar
#'
#' Bars of -log10(FDR), colored by direction (GSEA) or a single hue (ORA).
#'
#' @inheritParams fig_enrich_dot
#' @return a ggplot2 object
#' @export
fig_enrich_bar <- function(df, n = 20,
                           col_up = "#C0392B", col_down = "#2980B9",
                           mode = c("exploration", "publication")) {
  mode <- match.arg(mode)
  d <- enrich_plot_df(df, n)
  d$lbl <- clean_term(d$term)
  d$nlp <- -log10(d$padj + 1e-300)
  d <- d[order(d$nlp), ]
  d$lbl <- factor(d$lbl, levels = d$lbl)

  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$nlp, y = .data$lbl,
                                       fill = .data$direction)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_y_discrete(labels = function(x) wrap_label(x)) +
    ggplot2::labs(x = expression(-log[10] ~ FDR), y = NULL, fill = NULL) +
    fig_theme(mode)

  if (d$is_gsea[1]) {
    p <- p + ggplot2::scale_fill_manual(
      values = c("Up" = col_up, "Down" = col_down))
  } else {
    p <- p + ggplot2::scale_fill_manual(values = c("Enriched" = "#1D9E75")) +
      ggplot2::theme(legend.position = "none")
  }
  p
}

#' GSEA running-enrichment curve
#'
#' Wraps [fgsea::plotEnrichment()] for a single pathway, restyled to match
#' RNAflow's theme.
#'
#' @param res DE results data.frame
#' @param pathway_genes character vector of genes in the pathway
#' @param rank_by ranking metric passed to [rank_genes()]
#' @param title plot title (e.g. the pathway name)
#' @param line_color color of the running enrichment line
#' @param mode "exploration" or "publication"
#' @return a ggplot2 object
#' @export
fig_gsea_curve <- function(res, pathway_genes, rank_by = "stat",
                           title = NULL, line_color = "#1D9E75",
                           mode = c("exploration", "publication")) {
  mode <- match.arg(mode)
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    stop("Package 'fgsea' is required for the enrichment curve. ",
         "Install with: BiocManager::install('fgsea')", call. = FALSE)
  }
  ranks <- rank_genes(res, by = rank_by)
  genes <- intersect(as.character(pathway_genes), names(ranks))
  if (length(genes) < 2) {
    stop("Fewer than 2 pathway genes found in the ranked list.", call. = FALSE)
  }
  p <- fgsea::plotEnrichment(genes, ranks)
  # Restyle fgsea's default layers by geom type (robust across versions):
  # the running-ES line -> brand color, the gene-hit ticks -> grey.
  for (i in seq_along(p$layers)) {
    geom <- class(p$layers[[i]]$geom)[1]
    if (geom == "GeomLine") {
      p$layers[[i]]$aes_params$colour <- line_color
      p$layers[[i]]$aes_params$linewidth <- 0.8
    } else if (geom == "GeomSegment") {
      p$layers[[i]]$aes_params$colour <- "#BDC3C7"
    }
  }
  p +
    ggplot2::labs(title = clean_term(title %||% ""),
                  x = "Rank in ordered gene list",
                  y = "Enrichment score") +
    fig_theme(mode)
}
