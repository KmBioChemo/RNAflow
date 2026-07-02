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
    :root{--rf-accent:#1D9E75;--rf-accent-dark:#157a5b;--rf-ink:#1f2d3a;
      --rf-body:#46545f;--rf-muted:#7c8a94;--rf-border:#e6ebe9;--rf-surface-2:#f8faf9;}
    *{box-sizing:border-box;}
    body{font-family:'Inter',-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
      color:var(--rf-body);max-width:920px;margin:0 auto;padding:0 22px 56px;
      line-height:1.55;background:#fff;-webkit-font-smoothing:antialiased;}
    h1,h2,h3{color:var(--rf-ink);letter-spacing:-0.01em;}
    h2{font-size:1.28rem;margin-top:34px;padding-bottom:6px;
      border-bottom:1px solid var(--rf-border);}
    h3{font-size:1.02rem;color:#34495e;}
    a{color:var(--rf-accent-dark);text-decoration:none;}
    .report-header{margin:30px 0 8px;padding:22px 24px;border-radius:14px;color:#fff;
      background:linear-gradient(120deg,#1f2d3a 0%,#1D9E75 100%);
      box-shadow:0 6px 20px rgba(16,40,34,.12);}
    .report-header .rh-title{font-size:1.6rem;font-weight:800;letter-spacing:-0.02em;}
    .report-header .rh-sub{font-size:.92rem;opacity:.9;margin-top:4px;}
    .cards{display:flex;flex-wrap:wrap;gap:12px;margin:14px 0 6px;}
    .card{flex:1 1 150px;background:var(--rf-surface-2);border:1px solid var(--rf-border);
      border-radius:12px;padding:14px 16px;}
    .card .c-val{font-size:1.5rem;font-weight:750;color:var(--rf-ink);line-height:1.1;}
    .card .c-lbl{font-size:.72rem;font-weight:650;text-transform:uppercase;
      letter-spacing:.5px;color:var(--rf-muted);margin-top:4px;}
    table{border-collapse:collapse;width:100%;margin:10px 0;font-size:.86rem;}
    th,td{border:none;border-bottom:1px solid var(--rf-border);padding:7px 11px;text-align:left;}
    thead th{background:var(--rf-surface-2);color:var(--rf-ink);font-weight:650;
      border-bottom:2px solid #d6ddda;}
    tbody tr:nth-child(even){background:#fbfcfc;}
    pre{background:#0f172a;color:#e2e8f0;border:1px solid #0b1220;border-radius:10px;
      padding:14px 16px;overflow-x:auto;font-size:12.5px;line-height:1.5;
      font-family:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;}
    pre code{color:inherit;background:none;}
    .caption{color:var(--rf-muted);font-size:.8rem;margin:2px 0 16px;text-align:center;}
    .muted{color:var(--rf-muted);font-size:.85rem;}
    .badge{display:inline-block;background:#eafaf3;color:#0f7a54;border-radius:6px;
      padding:2px 9px;margin-right:6px;font-size:.8rem;font-weight:600;}
    .section{margin-top:8px;}
    .figure-block{margin:14px 0;}
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

  # ---- Overview cards (headline numbers) ----
  n_ann <- if (!is.null(meta)) max(0, ncol(meta) - 1) else 0
  gxs   <- if (!is.null(counts)) sprintf("%d &times; %d", nrow(counts), ncol(counts))
           else "&mdash;"
  o_card <- function(val, lbl) h$div(class = "card",
    h$div(class = "c-val", htmltools::HTML(val)), h$div(class = "c-lbl", lbl))
  overview_cards <- h$div(class = "cards",
    o_card(tools::toTitleCase(organism), "Organism"),
    o_card(gxs, "Genes &times; samples"),
    o_card(as.character(length(store)), "Contrasts"),
    o_card(as.character(n_ann), "Annotations"))

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
      h$div(class = "figure-block",
            h$h3(lbl, style = "font-size:15px;color:#34495e;"),
            embed_fig(p, 6, 4),
            h$div(class = "caption",
                  sprintf("Volcano plot -- %s (padj < 0.05, |log2FC| > 1).", lbl)))
    })
  }

  # ---- Cross-contrast comparison ----
  compare <- NULL
  if (length(store) >= 2) {
    contrasts <- contrast_store_results(store)
    hm <- tryCatch(fig_lfc_heatmap(contrasts, gene_src = "sig_union",
                                   n_genes = 40), error = function(e) NULL)
    if (!is.null(hm)) compare <- section(
      h$h2("Cross-contrast signature"),
      h$div(class = "figure-block", embed_fig(hm, 6, 6),
            h$div(class = "caption",
                  "log2 fold-change across contrasts for the union of ",
                  "significant genes (top 40).")))
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
    h$div(
      class = "report-header",
      h$div(class = "rh-title", title),
      h$div(class = "rh-sub",
            sprintf("Generated by RNAflow v%s%s", ver,
                    if (!is.null(generated)) paste0(" - ", generated) else ""))
    ),

    section(h$h2("Overview"), overview_cards, overview),
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
