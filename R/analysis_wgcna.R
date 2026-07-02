#' Weighted gene co-expression network analysis (WGCNA)
#'
#' Pure wrappers around \pkg{WGCNA} for the co-expression module: expression
#' matrix preparation, soft-threshold selection, blockwise module detection,
#' module-trait correlation, and intramodular hub genes. No Shiny dependency.
#'
#' Conventions: `counts_norm` is the normalized (e.g. vst) matrix with genes
#' in rows and samples in columns (as elsewhere in RNAflow). WGCNA wants the
#' transpose, so [wgcna_datexpr()] returns a samples x genes matrix
#' (`datExpr`) that the other functions consume.
#'
#' @name analysis_wgcna
#' @keywords internal
NULL

#' Run an expression while WGCNA's `cor` shadows `stats::cor`
#'
#' Several WGCNA routines call `do.call("cor", ...)` with WGCNA-specific
#' arguments (`weights.x`, `cosine`, ...). When the WGCNA package is not on
#' the search path (the normal case inside this package / a Shiny app), that
#' bare lookup resolves to `stats::cor` and errors. Binding `cor` to
#' `WGCNA::cor` on the global environment for the duration of the call is the
#' documented workaround; we restore the previous binding on exit.
#'
#' @param expr an expression to evaluate
#' @keywords internal
with_wgcna_cor <- function(expr) {
  ge  <- globalenv()
  had <- exists("cor", envir = ge, inherits = FALSE)
  old <- if (had) get("cor", envir = ge) else NULL
  assign("cor", WGCNA::cor, envir = ge)
  on.exit({
    if (had) assign("cor", old, envir = ge)
    else if (exists("cor", envir = ge, inherits = FALSE)) rm("cor", envir = ge)
  }, add = TRUE)
  force(expr)
}

#' Build the WGCNA expression matrix
#'
#' Selects the most variable genes and transposes to samples x genes.
#'
#' @param counts_norm normalized matrix (genes x samples)
#' @param n_genes number of top-variance genes to keep
#' @return a samples x genes numeric matrix (`datExpr`)
#' @export
wgcna_datexpr <- function(counts_norm, n_genes = 3000) {
  if (is.null(counts_norm) || !is.matrix(counts_norm) && !is.data.frame(counts_norm)) {
    stop("`counts_norm` must be a normalized genes x samples matrix.",
         call. = FALSE)
  }
  m <- as.matrix(counts_norm)
  if (ncol(m) < 4) {
    stop("WGCNA needs at least 4 samples; got ", ncol(m), ".", call. = FALSE)
  }
  v <- matrixStats::rowVars(m)
  keep <- order(v, decreasing = TRUE)[seq_len(min(n_genes, nrow(m)))]
  datExpr <- t(m[sort(keep), , drop = FALSE])

  # Flag/drop genes or samples with too many missing or zero-variance entries
  if (requireNamespace("WGCNA", quietly = TRUE)) {
    gsg <- with_wgcna_cor(WGCNA::goodSamplesGenes(datExpr, verbose = 0))
    if (!gsg$allOK) {
      nb_g <- sum(!gsg$goodGenes); nb_s <- sum(!gsg$goodSamples)
      message("goodSamplesGenes flagged ", nb_g, " gene(s) and ", nb_s,
              " sample(s) as low quality; removing them.")
      datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
    }
  }
  datExpr
}

#' Soft-threshold (power) selection
#'
#' @param datExpr samples x genes matrix from [wgcna_datexpr()]
#' @param powers candidate soft-thresholding powers
#' @param network_type "signed" (default), "unsigned", or "signed hybrid"
#' @param rsq_cut scale-free topology R^2 target used to suggest a power
#' @return a list with `fit_indices` (data.frame) and `suggested` (numeric)
#' @export
wgcna_pick_power <- function(datExpr, powers = 1:20,
                             network_type = "signed", rsq_cut = 0.8) {
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package 'WGCNA' is required. ",
         "Install with: BiocManager::install('WGCNA')", call. = FALSE)
  }
  sft <- with_wgcna_cor(WGCNA::pickSoftThreshold(
    datExpr, powerVector = powers, networkType = network_type, verbose = 0))
  fi <- sft$fitIndices
  # Suggest: first power reaching the R^2 target, else WGCNA's estimate.
  sft_rsq <- -sign(fi$slope) * fi$SFT.R.sq
  hit <- which(sft_rsq >= rsq_cut)
  fallback <- FALSE
  suggested <- if (length(hit) > 0) {
    fi$Power[hit[1]]
  } else if (!is.na(sft$powerEstimate)) {
    sft$powerEstimate
  } else {
    # Scale-free fit never reached the target (common with small sample sizes):
    # fall back to WGCNA's recommended default power for this sample size.
    fallback <- TRUE
    wgcna_default_power(nrow(datExpr), network_type)
  }
  list(fit_indices = fi, suggested = suggested, fallback = fallback,
       n_samples = nrow(datExpr), network_type = network_type, rsq_cut = rsq_cut)
}

#' WGCNA's recommended default soft-thresholding power
#'
#' Used when scale-free topology fit does not reach the target (e.g. small
#' sample sizes). Values follow the WGCNA FAQ recommendations.
#'
#' @param n_samples number of samples
#' @param network_type network type
#' @return a numeric power
#' @keywords internal
wgcna_default_power <- function(n_samples, network_type = "signed") {
  signed <- !identical(network_type, "unsigned")
  if (n_samples < 20)      if (signed) 18 else 9
  else if (n_samples < 30) if (signed) 16 else 8
  else if (n_samples < 40) if (signed) 14 else 7
  else                     if (signed) 12 else 6
}

#' Detect co-expression modules
#'
#' @param datExpr samples x genes matrix
#' @param power soft-thresholding power
#' @param network_type "signed" (default), "unsigned", "signed hybrid"
#' @param min_module_size minimum module size
#' @param merge_cut_height dendrogram cut for merging close modules
#' @param deep_split sensitivity of the tree cut (0-4)
#' @return a list: `modules` (named vector gene -> color), `MEs` (module
#'   eigengenes, samples x modules), `power`, `n_samples`, `datExpr`,
#'   `dendro` (the blockwise dendrogram)
#' @export
run_wgcna <- function(datExpr, power, network_type = "signed",
                      min_module_size = 30, merge_cut_height = 0.25,
                      deep_split = 2) {
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package 'WGCNA' is required. ",
         "Install with: BiocManager::install('WGCNA')", call. = FALSE)
  }
  # Keep the TOM type consistent with the chosen network type
  tom_type <- if (identical(network_type, "unsigned")) "unsigned" else "signed"
  net <- with_wgcna_cor(WGCNA::blockwiseModules(
    datExpr, power = power, networkType = network_type, TOMType = tom_type,
    minModuleSize = min_module_size, mergeCutHeight = merge_cut_height,
    deepSplit = deep_split, numericLabels = TRUE,
    maxBlockSize = ncol(datExpr) + 1, pamRespectsDendro = FALSE,
    saveTOMs = FALSE, randomSeed = 12345, verbose = 0))

  colors <- WGCNA::labels2colors(net$colors)
  names(colors) <- colnames(datExpr)
  MEs <- with_wgcna_cor(WGCNA::moduleEigengenes(datExpr, colors)$eigengenes)
  MEs <- WGCNA::orderMEs(MEs)
  rownames(MEs) <- rownames(datExpr)

  list(modules = colors, MEs = MEs, power = power,
       network_type = network_type, n_samples = nrow(datExpr),
       datExpr = datExpr,
       dendro = if (length(net$dendrograms)) net$dendrograms[[1]] else NULL,
       block_genes = net$blockGenes)
}

#' Build a numeric trait matrix from sample metadata
#'
#' Each annotation column is expanded into one indicator (0/1) column per
#' level, so module eigengenes can be correlated against every group.
#'
#' @param metadata data.frame (column 1 = sample ID, rest = annotations)
#' @param samples sample IDs to keep / order by (e.g. `rownames(datExpr)`)
#' @return a samples x traits numeric matrix
#' @export
build_traits <- function(metadata, samples) {
  if (is.null(metadata) || ncol(metadata) < 2) {
    stop("Metadata must have a sample column plus at least one annotation.",
         call. = FALSE)
  }
  samp_col <- colnames(metadata)[1]
  md <- metadata[match(samples, metadata[[samp_col]]), , drop = FALSE]
  ann <- colnames(metadata)[-1]
  cols <- list()
  for (a in ann) {
    f <- factor(as.character(md[[a]]))
    if (nlevels(f) < 2) next
    mm <- stats::model.matrix(~ 0 + f)
    colnames(mm) <- paste0(a, ": ", levels(f))
    for (j in seq_len(ncol(mm))) cols[[colnames(mm)[j]]] <- mm[, j]
  }
  if (length(cols) == 0) {
    stop("No annotation column had >= 2 levels to build traits from.",
         call. = FALSE)
  }
  out <- do.call(cbind, cols)
  rownames(out) <- samples
  out
}

#' Module-trait correlation
#'
#' @param MEs module eigengenes (samples x modules) from [run_wgcna()]
#' @param traits numeric trait matrix from [build_traits()]
#' @return a list: `cor` (modules x traits), `p` (raw p-values), `padj`
#'   (Benjamini-Hochberg across the whole matrix), `n` (samples)
#' @details Many correlations are tested at once, so `padj` applies BH
#'   correction across all module x trait cells; prefer it over raw `p`.
#' @export
module_trait_cor <- function(MEs, traits) {
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package 'WGCNA' is required. ",
         "Install with: BiocManager::install('WGCNA')", call. = FALSE)
  }
  common <- intersect(rownames(MEs), rownames(traits))
  if (length(common) < 3) {
    stop("Fewer than 3 samples shared between eigengenes and traits.",
         call. = FALSE)
  }
  MEs <- MEs[common, , drop = FALSE]
  traits <- traits[common, , drop = FALSE]
  cmat <- stats::cor(MEs, traits, use = "pairwise.complete.obs")
  pmat <- WGCNA::corPvalueStudent(cmat, length(common))
  padj <- matrix(stats::p.adjust(as.vector(pmat), method = "BH"),
                 nrow = nrow(pmat), dimnames = dimnames(pmat))
  list(cor = cmat, p = pmat, padj = padj, n = length(common))
}

#' Intramodular hub genes
#'
#' Ranks the genes of a module by their module membership (signed kME).
#'
#' @param wg the list returned by [run_wgcna()]
#' @param module module color
#' @param n number of hub genes to return
#' @return a data.frame (gene, kME) sorted by descending kME
#' @export
hub_genes <- function(wg, module, n = 20) {
  if (!requireNamespace("WGCNA", quietly = TRUE)) {
    stop("Package 'WGCNA' is required. ",
         "Install with: BiocManager::install('WGCNA')", call. = FALSE)
  }
  if (!module %in% wg$modules) {
    stop("Module '", module, "' not found.", call. = FALSE)
  }
  genes <- names(wg$modules)[wg$modules == module]
  kme <- with_wgcna_cor(WGCNA::signedKME(wg$datExpr, wg$MEs))
  col <- paste0("kME", module)
  if (!col %in% colnames(kme)) {
    stop("kME column for module '", module, "' not found.", call. = FALSE)
  }
  sc <- kme[genes, col]
  ord <- order(sc, decreasing = TRUE)
  data.frame(gene = genes[ord], kME = round(sc[ord], 4),
             row.names = NULL, stringsAsFactors = FALSE)[seq_len(min(n, length(genes))), ]
}

#' Module gene lists
#'
#' @param wg the list returned by [run_wgcna()]
#' @param exclude_grey drop the unassigned "grey" module (default TRUE)
#' @return a named list of character vectors (module color -> genes), ready
#'   to feed [run_ora()] for module-to-pathway enrichment
#' @export
module_gene_list <- function(wg, exclude_grey = TRUE) {
  gl <- split(names(wg$modules), wg$modules)
  if (isTRUE(exclude_grey)) gl[["grey"]] <- NULL
  gl
}

#' Per-module summary table
#'
#' @param wg the list returned by [run_wgcna()]
#' @return a data.frame (module, n_genes) sorted by size, grey last
#' @export
module_summary <- function(wg) {
  tab <- table(wg$modules)
  df <- data.frame(module = names(tab), n_genes = as.integer(tab),
                   stringsAsFactors = FALSE)
  is_grey <- df$module == "grey"
  df <- rbind(df[!is_grey, ][order(-df$n_genes[!is_grey]), ], df[is_grey, ])
  rownames(df) <- NULL
  df
}
