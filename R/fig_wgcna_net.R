#' WGCNA module co-expression network figure
#'
#' Draws a module as a co-expression network: the top genes by module
#' membership (kME) as nodes, edges where the pairwise correlation exceeds a
#' threshold. Hub genes sit at the center. Built on \pkg{ggraph}.
#'
#' @name fig_wgcna_net
#' @keywords internal
NULL

#' Module co-expression network
#'
#' @param wg the list returned by [run_wgcna()]
#' @param module module color
#' @param n number of top-kME genes to include
#' @param min_cor minimum absolute correlation to draw an edge
#' @param label_n number of top hub genes to label
#' @param mode "exploration" or "publication"
#' @return a ggplot2 / ggraph object
#' @export
fig_module_network <- function(wg, module, n = 30, min_cor = 0.4,
                               label_n = 12,
                               mode = c("exploration", "publication")) {
  mode <- match.arg(mode)
  for (p in c("ggraph", "igraph", "tidygraph")) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop("Package '", p, "' is required for the module network. ",
           "Install with: install.packages('", p, "')", call. = FALSE)
    }
  }
  hubs <- hub_genes(wg, module, n = n)
  genes <- hubs$gene
  if (length(genes) < 3) {
    stop("Fewer than 3 genes in module '", module, "' to build a network.",
         call. = FALSE)
  }
  expr <- wg$datExpr[, genes, drop = FALSE]
  cmat <- stats::cor(expr)

  # Upper-triangle edges above threshold
  k <- length(genes)
  ij <- which(upper.tri(cmat) & abs(cmat) >= min_cor, arr.ind = TRUE)
  edges <- if (nrow(ij) > 0)
    data.frame(from = ij[, 1], to = ij[, 2],
               weight = abs(cmat[ij]), sign = sign(cmat[ij]))
  else data.frame(from = integer(0), to = integer(0),
                  weight = numeric(0), sign = numeric(0))

  ord <- order(hubs$kME, decreasing = TRUE)
  nodes <- data.frame(name = genes, kME = hubs$kME, stringsAsFactors = FALSE)
  nodes$lbl <- NA_character_
  nodes$lbl[ord[seq_len(min(label_n, k))]] <- genes[ord[seq_len(min(label_n, k))]]

  g <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = FALSE)

  p <- ggraph::ggraph(g, layout = "fr")
  if (nrow(edges) > 0) {
    p <- p + ggraph::geom_edge_link(
      ggplot2::aes(edge_alpha = .data$weight),
      edge_colour = "#B7C0C9", edge_width = 0.4, show.legend = FALSE)
  }
  p +
    ggraph::geom_node_point(
      ggplot2::aes(size = .data$kME, fill = .data$kME),
      shape = 21, colour = "white", stroke = 0.4) +
    ggraph::geom_node_text(
      ggplot2::aes(label = .data$lbl), size = 2.6, repel = TRUE,
      fontface = "italic", family = "", colour = "#2C3E50",
      na.rm = TRUE, max.overlaps = 30,
      box.padding = 0.5, point.padding = 0.3, min.segment.length = 0.4,
      segment.colour = "#C0C7CE", segment.size = 0.25,
      bg.colour = "white", bg.r = 0.14) +
    scale_fill_omics_c("batlow") +
    ggplot2::scale_size_continuous(range = c(2.5, 10), guide = "none") +
    ggplot2::labs(fill = "kME",
                  title = sprintf("Module %s co-expression network", module)) +
    ggraph::theme_graph(base_family = "") +
    ggplot2::theme(legend.position = "right",
                   plot.title = ggplot2::element_text(face = "bold", size = 12),
                   plot.background = ggplot2::element_rect(fill = "white", colour = NA))
}
