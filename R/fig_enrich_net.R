#' Network and distribution enrichment figures
#'
#' Additional enrichment visuals: an enrichment map (pathways as a network
#' linked by shared genes) and a GSEA ridgeline (per-pathway distributions of
#' the ranking metric). Built on \pkg{ggraph} / \pkg{igraph} and
#' \pkg{ggridges}. Publication-oriented aesthetics via the shared theme.
#'
#' @name fig_enrich_net
#' @keywords internal
NULL

#' Per-pathway gene lists from a GSEA or ORA table
#'
#' @param df a data.frame from [run_gsea()] (uses `leadingEdge`) or
#'   [run_ora()] (uses `geneID`)
#' @param n keep the top `n` terms by padj
#' @return a list with `nodes` (term, score, size, is_gsea) and `genes`
#'   (named list of character vectors)
#' @keywords internal
enrich_network_data <- function(df, n = 30) {
  if (!is.data.frame(df) || nrow(df) == 0) {
    stop("No enrichment results to plot.", call. = FALSE)
  }
  df <- df[order(df$padj), , drop = FALSE]
  df <- utils::head(df, max(2L, as.integer(n)))
  is_gsea <- "NES" %in% colnames(df)
  if (is_gsea) {
    term  <- as.character(df$pathway)
    score <- as.numeric(df$NES)
    size  <- as.numeric(df$size %||% NA)
    genes <- if (is.list(df$leadingEdge)) lapply(df$leadingEdge, as.character)
             else strsplit(as.character(df$leading_edge %||% ""), ",\\s*")
  } else {
    term  <- as.character(df$Description)
    score <- -log10(as.numeric(df$padj) + 1e-300)
    size  <- as.numeric(df$Count %||% NA)
    genes <- strsplit(as.character(df$geneID %||% ""), "/")
  }
  names(genes) <- term
  list(
    nodes = data.frame(term = term, score = score, size = size,
                       is_gsea = is_gsea, stringsAsFactors = FALSE),
    genes = lapply(genes, function(g) unique(g[nzchar(g)]))
  )
}

#' Enrichment map (pathway network)
#'
#' Nodes are enriched terms, edges connect terms sharing genes (Jaccard
#' similarity). Node color encodes NES (GSEA) or -log10(FDR) (ORA), node size
#' encodes gene-set size. A compact overview of how the enriched
#' biology clusters.
#'
#' @param df a data.frame from [run_gsea()] or [run_ora()]
#' @param n number of top terms (by padj) to include
#' @param min_similarity minimum Jaccard similarity to draw an edge
#' @param label_n number of top nodes (by score magnitude) to label
#' @param mode "exploration" or "publication"
#' @return a ggplot2 / ggraph object
#' @export
fig_enrich_map <- function(df, n = 30, min_similarity = 0.2, label_n = 12,
                           mode = c("exploration", "publication")) {
  mode <- match.arg(mode)
  for (p in c("ggraph", "igraph", "tidygraph")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Package '", p, "' is required for the enrichment map. ",
           "Install with: install.packages('", p, "')", call. = FALSE)
    }
  }
  nd <- enrich_network_data(df, n)
  nodes <- nd$nodes
  nodes$label <- clean_term(nodes$term)
  genes <- nd$genes

  # Jaccard edges
  k <- nrow(nodes)
  edges <- data.frame(from = integer(0), to = integer(0), weight = numeric(0))
  if (k >= 2) {
    for (i in seq_len(k - 1)) for (j in (i + 1):k) {
      a <- genes[[i]]; b <- genes[[j]]
      u <- length(union(a, b))
      w <- if (u > 0) length(intersect(a, b)) / u else 0
      if (w >= min_similarity) edges <- rbind(edges,
        data.frame(from = i, to = j, weight = w))
    }
  }

  # Label only the strongest nodes to avoid clutter
  ord <- order(abs(nodes$score), decreasing = TRUE)
  nodes$show <- FALSE
  nodes$show[utils::head(ord, max(0L, as.integer(label_n)))] <- TRUE
  nodes$lbl <- ifelse(nodes$show, nodes$label, NA_character_)

  g <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = FALSE)
  is_gsea <- nodes$is_gsea[1]
  score_lab <- if (is_gsea) "NES" else expression(-log[10] ~ FDR)

  p <- ggraph::ggraph(g, layout = "fr")
  if (nrow(edges) > 0) {
    p <- p + ggraph::geom_edge_link(
      ggplot2::aes(edge_alpha = .data$weight),
      edge_colour = "#9AA6B2", edge_width = 0.5, show.legend = FALSE)
  }
  fill_scale <- if (is_gsea) {
    lim <- max(abs(nodes$score), na.rm = TRUE)
    lim <- if (is.finite(lim) && lim > 0) lim else 1
    scale_fill_omics_div("vik", midpoint = 0, limits = c(-lim, lim))
  } else scale_fill_omics_c("batlow")

  p +
    ggraph::geom_node_point(
      ggplot2::aes(size = .data$size, fill = .data$score),
      shape = 21, colour = "white", stroke = 0.4) +
    ggraph::geom_node_text(
      ggplot2::aes(label = .data$lbl), size = 2.6, repel = TRUE,
      family = "", colour = "#2C3E50", na.rm = TRUE, max.overlaps = 30,
      box.padding = 0.6, point.padding = 0.4, min.segment.length = 0.4,
      segment.colour = "#C0C7CE", segment.size = 0.25,
      bg.colour = "white", bg.r = 0.14) +
    fill_scale +
    ggplot2::scale_size_continuous(range = c(3, 11), guide = "none") +
    ggplot2::labs(fill = score_lab) +
    ggraph::theme_graph(base_family = "") +
    ggplot2::theme(legend.position = "right",
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA))
}

#' GSEA ridgeline plot
#'
#' For the top pathways, the distribution of the gene-level ranking metric
#' across each pathway's members, drawn as stacked density ridges colored by
#' NES. Shows *how* each set is shifted in the ranked list, not just a single
#' score.
#'
#' @param res DE results data.frame
#' @param gene_sets named list of gene sets (from [get_gene_sets()])
#' @param gsea a data.frame from [run_gsea()] (for ordering / NES coloring)
#' @param n number of top pathways (by padj) to show
#' @param rank_by ranking metric passed to [rank_genes()]
#' @param mode "exploration" or "publication"
#' @return a ggplot2 object
#' @export
fig_gsea_ridge <- function(res, gene_sets, gsea, n = 15, rank_by = "stat",
                           mode = c("exploration", "publication")) {
  mode <- match.arg(mode)
  if (!requireNamespace("ggridges", quietly = TRUE)) {
    stop("Package 'ggridges' is required for the ridgeline plot. ",
         "Install with: install.packages('ggridges')", call. = FALSE)
  }
  if (!is.data.frame(gsea) || nrow(gsea) == 0) {
    stop("No GSEA results to plot.", call. = FALSE)
  }
  ranks <- rank_genes(res, by = rank_by)
  top <- utils::head(gsea[order(gsea$padj), , drop = FALSE], max(2L, as.integer(n)))

  parts <- lapply(seq_len(nrow(top)), function(i) {
    pw <- top$pathway[i]
    genes <- intersect(gene_sets[[pw]], names(ranks))
    if (length(genes) < 3) return(NULL)
    data.frame(pathway = pw, value = unname(ranks[genes]),
               NES = top$NES[i], stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, parts)
  if (is.null(df) || nrow(df) == 0) {
    stop("Not enough overlapping genes to build ridges.", call. = FALSE)
  }
  # Order pathways by NES (bottom = most negative)
  lev <- unique(df[order(df$NES), "pathway"])
  df$lbl <- factor(clean_term(df$pathway),
                   levels = clean_term(lev))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$value, y = .data$lbl,
                                   fill = .data$NES)) +
    ggridges::geom_density_ridges(scale = 2.2, rel_min_height = 0.01,
                                  colour = "white", linewidth = 0.3,
                                  alpha = 0.95) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "#7F8C8D", linewidth = 0.4) +
    scale_fill_omics_div("vik", midpoint = 0,
                         limits = {
                           lim <- max(abs(df$NES), na.rm = TRUE)
                           lim <- if (is.finite(lim) && lim > 0) lim else 1
                           c(-lim, lim)
                         }) +
    ggplot2::labs(x = "Gene ranking metric", y = NULL, fill = "NES") +
    fig_theme(mode)
}
