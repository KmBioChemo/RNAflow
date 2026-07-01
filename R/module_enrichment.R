#' Module-level functional enrichment
#'
#' Runs over-representation analysis on each WGCNA module and visualizes the
#' result as a modules x pathways dotplot (compareCluster-style), so you can
#' see at a glance which biology each co-expression module captures. Reuses
#' the phase-3 enrichment layer.
#'
#' @name module_enrichment
#' @keywords internal
NULL

#' Enrich every module against a pathway database
#'
#' @param wg the list returned by [run_wgcna()]
#' @param organism one of "human", "mouse", "rat"
#' @param db "GO", "KEGG", or "Reactome"
#' @param ont GO ontology when `db == "GO"`
#' @param n_per number of top terms (by padj) to keep per module
#' @param padj_cutoff adjusted p-value cutoff passed to the ORA
#' @param min_size,max_size gene-set size filters
#' @return a tidy data.frame combining the per-module ORA tables, with an
#'   extra `module` column (only modules with >= 1 enriched term are kept)
#' @export
enrich_modules <- function(wg, organism, db = c("GO", "KEGG", "Reactome"),
                           ont = c("BP", "MF", "CC"), n_per = 5,
                           padj_cutoff = 0.1, min_size = 10, max_size = 500) {
  db  <- match.arg(db)
  ont <- match.arg(ont)
  gl <- module_gene_list(wg, exclude_grey = TRUE)
  if (length(gl) == 0) stop("No non-grey modules to enrich.", call. = FALSE)
  universe <- colnames(wg$datExpr)

  parts <- list()
  for (m in names(gl)) {
    er <- tryCatch(
      run_ora(gl[[m]], organism, db = db,
              ont = if (db == "GO") ont else "BP", universe = universe,
              padj_cutoff = padj_cutoff, min_size = min_size,
              max_size = max_size),
      error = function(e) NULL)
    if (is.null(er) || nrow(er) == 0) next
    er <- utils::head(er[order(er$padj), , drop = FALSE], max(1L, as.integer(n_per)))
    er$module <- m
    parts[[m]] <- er
  }
  if (length(parts) == 0) {
    stop("No enriched terms found in any module. Try a looser cutoff or a ",
         "different database.", call. = FALSE)
  }
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

#' Module x pathway enrichment dotplot
#'
#' @param combined a data.frame from [enrich_modules()]
#' @param max_terms maximum number of distinct terms (rows) to display
#' @param mode "exploration" or "publication"
#' @return a ggplot2 object
#' @export
fig_module_enrichment <- function(combined, max_terms = 25,
                                  mode = c("exploration", "publication")) {
  mode <- match.arg(mode)
  if (!is.data.frame(combined) || nrow(combined) == 0 ||
      !all(c("module", "Description", "Count", "padj") %in% colnames(combined))) {
    stop("`combined` must come from enrich_modules() (non-empty).", call. = FALSE)
  }
  df <- combined
  df$nlp <- -log10(as.numeric(df$padj) + 1e-300)
  df$Count <- as.numeric(df$Count)

  # Keep the strongest terms overall, then order rows by peak module
  best <- tapply(df$nlp, df$Description, max)
  keep_terms <- names(sort(best, decreasing = TRUE))[
    seq_len(min(max_terms, length(best)))]
  df <- df[df$Description %in% keep_terms, , drop = FALSE]

  mod_levels <- unique(df$module)
  # order terms so related modules cluster: by (peak module index, -nlp)
  peak_mod <- vapply(split(df, df$Description),
                     function(d) match(d$module[which.max(d$nlp)], mod_levels),
                     numeric(1))
  term_order <- names(sort(peak_mod, decreasing = TRUE))
  df$Description <- factor(df$Description, levels = term_order)
  df$module <- factor(df$module, levels = mod_levels)

  # Color module axis labels by their WGCNA color when ggtext is available
  use_md <- requireNamespace("ggtext", quietly = TRUE)
  x_labels <- if (use_md) {
    stats::setNames(sprintf("<span style='color:%s'>**%s**</span>",
                            mod_levels, mod_levels), mod_levels)
  } else ggplot2::waiver()

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$module, y = .data$Description,
                                        size = .data$Count, colour = .data$nlp)) +
    ggplot2::geom_point() +
    scale_color_omics_c("batlow") +
    ggplot2::scale_size_continuous(range = c(2, 8), name = "Gene count") +
    ggplot2::scale_y_discrete(labels = function(x) wrap_label(x, 40)) +
    ggplot2::labs(x = NULL, y = NULL, colour = expression(-log[10] ~ FDR)) +
    fig_theme(mode)

  if (use_md) {
    p <- p + ggplot2::scale_x_discrete(labels = x_labels) +
      ggplot2::theme(
        axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1))
  } else {
    p <- p + ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
  p
}
