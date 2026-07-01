#' Methods paragraph generator
#'
#' Turns an analysis session into a short, citable Methods paragraph naming
#' the tools, their versions, and the parameters actually used -- a prose
#' companion to [generate_r_script()], ready to adapt for a manuscript.
#'
#' @name methods_text
#' @keywords internal
NULL

# " (v1.2.3)" if the package is installed, else "".
pkg_ver <- function(p) {
  v <- tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA)
  if (is.na(v)) "" else sprintf(" (v%s)", v)
}

# Human-readable label for a GSEA ranking metric.
rank_label <- function(x) {
  switch(x %||% "stat",
         stat = "the Wald statistic",
         signed_p = "the signed -log10 p-value",
         log2fc = "the log2 fold change",
         x)
}

#' Generate a Methods paragraph for an analysis
#'
#' @param project a project list (organism, contrasts, enrichment, wgcna)
#' @return a single character string of prose
#' @export
generate_methods_text <- function(project) {
  store <- project$contrasts %||% list()
  en <- project$enrichment
  wg <- project$wgcna
  organism <- project$organism %||% "the study organism"
  parts <- character(0)

  ## ---- Differential expression ----
  if (length(store) > 0) {
    p1 <- store[[1]]$params %||% list()
    covs <- p1$covariates
    min_count <- p1$min_count %||% 10
    alpha <- p1$alpha %||% 0.05
    shrink_used <- p1$shrink_used %||% (if (isTRUE(p1$shrink)) "apeglm" else "none")

    labels <- vapply(names(store), function(nm) {
      pp <- store[[nm]]$params %||% list()
      if (!is.null(pp$treated))
        sprintf("%s vs %s (variable '%s')", pp$treated, pp$reference, pp$design_var)
      else nm
    }, character(1))
    contrast_sentence <- if (length(labels) == 1)
      sprintf("The contrast tested was %s.", labels[1])
    else sprintf("The %d contrasts tested were: %s.", length(labels),
                 paste(labels, collapse = "; "))

    adjust <- if (length(covs))
      sprintf(", adjusting for %s", paste(covs, collapse = " and ")) else ""
    shrink_sentence <- if (shrink_used == "none")
      "Effect sizes were reported without shrinkage."
    else sprintf(paste0("Log2 fold-change estimates were shrunk with the '%s' ",
                        "estimator for ranking and visualization, while ",
                        "statistical inference (Wald test, ",
                        "Benjamini-Hochberg-adjusted) used the unshrunken model."),
                 shrink_used)

    parts <- c(parts, sprintf(
      paste0("Differential expression analysis was performed with DESeq2%s. ",
             "Genes with a total count below %s were removed before model ",
             "fitting%s. %s %s"),
      pkg_ver("DESeq2"), min_count, adjust, contrast_sentence, shrink_sentence))
  }

  ## ---- GSEA / ORA ----
  if (is.list(en) && length(en) && identical(en$method, "gsea")) {
    sub <- if (!is.null(en$subcollection))
      sprintf(" (subcollection %s)", en$subcollection) else ""
    parts <- c(parts, sprintf(
      paste0("Gene set enrichment analysis was performed with fgsea%s against ",
             "the MSigDB %s collection%s (retrieved with msigdbr%s), ranking ",
             "genes by %s."),
      pkg_ver("fgsea"), en$collection, sub, pkg_ver("msigdbr"),
      rank_label(en$rank_by)))
  } else if (is.list(en) && length(en) && identical(en$method, "ora")) {
    ont <- if (identical(en$db, "GO")) sprintf(" (%s ontology)", en$ont) else ""
    parts <- c(parts, sprintf(
      paste0("Over-representation analysis was performed with clusterProfiler%s ",
             "against %s%s, testing genes with adjusted p < %s and ",
             "|log2 fold change| > %s against the background of all tested genes."),
      pkg_ver("clusterProfiler"), en$db, ont, en$padj %||% 0.05, en$lfc %||% 1))
  }

  ## ---- WGCNA ----
  if (is.list(wg) && length(wg) && !is.null(wg$power)) {
    parts <- c(parts, sprintf(
      paste0("Weighted gene co-expression network analysis was performed with ",
             "WGCNA%s on the %s most variable genes, using a %s network with a ",
             "soft-thresholding power of %s; modules were detected with a ",
             "minimum size of %s and merged at a dissimilarity threshold of %s."),
      pkg_ver("WGCNA"), wg$n_genes %||% 3000, wg$network_type %||% "signed",
      wg$power, wg$min_module_size %||% 30, wg$merge_cut_height %||% 0.25))
  }

  if (length(parts) == 0) {
    parts <- sprintf(paste0("No analysis has been run yet. Once a contrast is ",
                            "computed for %s, this paragraph will summarize the ",
                            "methods used."), organism)
  }

  footer <- sprintf(
    "Analyses were carried out in R (v%s) with RNAflow%s.",
    getRversion(), pkg_ver("RNAflow"))
  paste(c(parts, footer), collapse = " ")
}
