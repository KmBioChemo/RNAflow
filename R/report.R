#' Standalone HTML analysis report
#'
#' Builds a self-contained HTML report of an RNAflow session using
#' \pkg{htmltools} -- figures are embedded as base64 data URIs, so the file
#' needs no external assets and no pandoc / Quarto toolchain. Fast, dependable
#' figures (DE volcanoes, cross-contrast comparison) are rendered inline; the
#' heavier enrichment / network steps are documented via the reproducible
#' script (see [generate_r_script()]).
#'
#' @name report
#' @keywords internal
NULL

#' Key package versions used in the session
#'
#' @return a data.frame (package, version) for packages that are installed
#' @keywords internal
session_manifest <- function() {
  pkgs <- c("RNAflow", "DESeq2", "SummarizedExperiment", "ggplot2", "fgsea",
            "clusterProfiler", "msigdbr", "ReactomePA", "WGCNA", "eulerr",
            "ComplexHeatmap", "pheatmap", "plotly", "crosstalk", "httr2",
            "decoupleR", "OmnipathR")
  rows <- lapply(pkgs, function(p) {
    v <- tryCatch(as.character(utils::packageVersion(p)),
                  error = function(e) NA_character_)
    data.frame(package = p, version = v, stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  df[!is.na(df$version), , drop = FALSE]
}

#' Embed any RNAflow figure as an <img> data URI
#' @keywords internal
embed_fig <- function(obj, w = 6, h = 4, dpi = 110) {
  tmp <- tempfile(fileext = ".png")
  on.exit(unlink(tmp), add = TRUE)
  ok <- tryCatch({ save_compare(obj, tmp, "png", w, h, dpi); TRUE },
                 error = function(e) FALSE)
  if (!ok || !file.exists(tmp)) {
    return(htmltools::tags$p(htmltools::tags$em("(figure could not be rendered)")))
  }
  htmltools::tags$img(src = xfun::base64_uri(tmp),
                      style = "max-width:100%;height:auto;margin:8px 0;")
}

#' Per-contrast DE summary counts
#' @keywords internal
de_counts <- function(res, padj_thr = 0.05, lfc_thr = 1) {
  ok <- !is.na(res$padj) & !is.na(res$log2FoldChange)
  up <- sum(ok & res$padj < padj_thr & res$log2FoldChange >  lfc_thr)
  dn <- sum(ok & res$padj < padj_thr & res$log2FoldChange < -lfc_thr)
  c(up = up, down = dn, total = nrow(res))
}

#' Build a standalone HTML report for an analysis session
#'
#' @param project a project list (organism, counts, metadata, contrasts store)
#' @param file output HTML path
#' @param title report title
#' @param generated optional timestamp string for the header
#' @return invisibly, the output path
#' @export
build_report_html <- function(project, file, title = "RNAflow analysis report",
                              generated = NULL) {
  organism <- project$organism %||% "unknown"
  store    <- project$contrasts %||% list()
  counts   <- project$counts
  meta     <- project$metadata
  ver <- tryCatch(as.character(utils::packageVersion("RNAflow")),
                  error = function(e) "dev")

  h <- htmltools::tags
  section <- function(...) h$div(class = "section", ...)

  css <- "
    body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
      color:#2C3E50;max-width:900px;margin:24px auto;padding:0 18px;line-height:1.5;}
    h1{color:#1D9E75;border-bottom:2px solid #1D9E75;padding-bottom:6px;}
    h2{color:#2C3E50;margin-top:28px;border-bottom:1px solid #e9ecef;padding-bottom:4px;}
    table{border-collapse:collapse;width:100%;margin:8px 0;font-size:14px;}
    th,td{border:1px solid #dee2e6;padding:5px 9px;text-align:left;}
    th{background:#f2f7f5;}
    pre{background:#f6f8fa;border:1px solid #e1e4e8;border-radius:6px;padding:12px;
      overflow-x:auto;font-size:12.5px;}
    .muted{color:#7F8C8D;font-size:13px;}
    .badge{display:inline-block;background:#eafaf3;color:#0f7a54;border-radius:4px;
      padding:1px 7px;margin-right:6px;font-size:13px;}
  "

  # ---- Overview table ----
  overview <- h$table(
    h$tr(h$th("Organism"), h$td(organism)),
    h$tr(h$th("Genes x samples"),
         h$td(if (!is.null(counts)) sprintf("%d x %d", nrow(counts), ncol(counts))
              else "not stored")),
    h$tr(h$th("Metadata annotations"),
         h$td(if (!is.null(meta)) paste(colnames(meta)[-1], collapse = ", ")
              else "not stored")),
    h$tr(h$th("Saved contrasts"), h$td(length(store)))
  )

  # ---- Contrasts summary table ----
  contrast_tbl <- NULL
  if (length(store) > 0) {
    rows <- lapply(names(store), function(lbl) {
      cnt <- de_counts(store[[lbl]]$results)
      h$tr(h$td(lbl), h$td(cnt["up"]), h$td(cnt["down"]), h$td(cnt["total"]))
    })
    contrast_tbl <- h$table(
      h$tr(h$th("Contrast"), h$th("Up"), h$th("Down"), h$th("Genes")), rows)
  }

  # ---- Per-contrast volcanoes ----
  volcanoes <- NULL
  if (length(store) > 0) {
    volcanoes <- lapply(names(store), function(lbl) {
      p <- tryCatch(fig_volcano(store[[lbl]]$results, lfc_thr = 1, padj_thr = 0.05,
                                n_label = 15, title = lbl),
                    error = function(e) NULL)
      if (is.null(p)) return(NULL)
      h$div(h$h3(lbl, style = "font-size:15px;color:#34495e;"),
            embed_fig(p, 6, 4))
    })
  }

  # ---- Cross-contrast comparison ----
  compare <- NULL
  if (length(store) >= 2) {
    contrasts <- contrast_store_results(store)
    hm <- tryCatch(fig_lfc_heatmap(contrasts, gene_src = "sig_union",
                                   n_genes = 40), error = function(e) NULL)
    if (!is.null(hm)) compare <- section(h$h2("Cross-contrast signature"),
                                         embed_fig(hm, 6, 6))
  }

  # ---- AI interpretation (optional) ----
  ai <- project$ai_interpretation
  ai_section <- NULL
  if (is.list(ai) && !is.null(ai$text) && nzchar(ai$text)) {
    ai_html <- tryCatch(
      htmltools::HTML(as.character(shiny::markdown(ai$text))),
      error = function(e) h$pre(h$code(ai$text)))
    prov <- sprintf("Generated with %s via the Anthropic Claude API%s%s.",
                    ai$model %||% "an LLM",
                    if (!is.null(ai$generated)) paste0(" on ", ai$generated) else "",
                    if (!is.null(ai$input_tokens) && !is.na(ai$input_tokens))
                      sprintf(" (%s in / %s out tokens%s)",
                              ai$input_tokens, ai$output_tokens,
                              if (!is.null(ai$cost_usd) && !is.na(ai$cost_usd))
                                sprintf(", ~$%.4f", ai$cost_usd) else "")
                    else "")
    ai_section <- section(
      h$h2("AI interpretation"),
      h$p(class = "muted",
          prov, " AI-generated text can be wrong -- treat it as a ",
          "hypothesis-generating draft, not a conclusion."),
      h$div(style = "background:#f8fbfa;border:1px solid #e1efe9;border-radius:6px;padding:2px 16px;",
            ai_html))
  }

  # ---- Reproducible script ----
  script <- generate_r_script(project, generated = generated)

  # ---- Session manifest ----
  man <- session_manifest()
  man_tbl <- h$table(
    h$tr(h$th("Package"), h$th("Version")),
    lapply(seq_len(nrow(man)),
           function(i) h$tr(h$td(man$package[i]), h$td(man$version[i]))))

  doc <- htmltools::tagList(
    h$head(h$style(htmltools::HTML(css)), h$title(title)),
    h$h1(title),
    h$p(class = "muted",
        sprintf("Generated by RNAflow v%s%s", ver,
                if (!is.null(generated)) paste0(" - ", generated) else "")),

    section(h$h2("Overview"), overview),
    if (!is.null(contrast_tbl))
      section(h$h2("Differential expression"),
              h$p(class = "muted", "Thresholds: padj < 0.05, |log2FC| > 1."),
              contrast_tbl, volcanoes),
    compare,
    ai_section,
    section(h$h2("Reproducible R script"),
            h$p(class = "muted",
                "A script that reproduces the pipeline with RNAflow (with the ",
                "package installed). DE calls reflect each saved contrast; ",
                "enrichment and network steps use the settings recorded during ",
                "the session when available, and documented example defaults ",
                "otherwise."),
            h$pre(h$code(script))),
    section(h$h2("Session"),
            h$p(h$span(class = "badge", paste0("RNAflow ", ver)),
                h$span(class = "badge", paste0("R ", getRversion()))),
            man_tbl,
            h$h3("sessionInfo()", style = "font-size:15px;margin-top:14px;"),
            h$pre(h$code(paste(utils::capture.output(utils::sessionInfo()),
                               collapse = "\n"))))
  )

  htmltools::save_html(doc, file = file, background = "white")
  invisible(file)
}
