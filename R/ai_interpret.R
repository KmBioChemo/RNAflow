#' AI-assisted biological interpretation
#'
#' Turns a differential-expression contrast (plus optional functional
#' enrichment) into a compact prompt and asks Anthropic's Claude API to write
#' a biological narrative. The prompt-building functions are pure and testable
#' without any network access; [call_claude()] is the only function that
#' touches the API (guarded by \pkg{httr2}).
#'
#' Only a compact summary leaves the machine -- top gene *names* with their
#' fold-changes and FDRs, and the top enrichment terms. Count matrices and
#' sample metadata are never sent.
#'
#' @name ai_interpret
#' @keywords internal
NULL

# Models offered in the UI (label -> API id). Opus 4.8 is the default.
AI_MODELS <- c(
  "Claude Opus 4.8 (best)"      = "claude-opus-4-8",
  "Claude Sonnet 5"             = "claude-sonnet-5",
  "Claude Haiku 4.5 (cheapest)" = "claude-haiku-4-5"
)

# USD per 1M tokens (input, output) for cost estimation.
MODEL_PRICING <- list(
  "claude-opus-4-8"   = c(in_ = 5,  out = 25),
  "claude-opus-4-7"   = c(in_ = 5,  out = 25),
  "claude-sonnet-5"   = c(in_ = 3,  out = 15),
  "claude-sonnet-4-6" = c(in_ = 3,  out = 15),
  "claude-haiku-4-5"  = c(in_ = 1,  out = 5)
)

# System prompt: keep the model rigorous and grounded.
RNAFLOW_AI_SYSTEM <- paste0(
  "You are a computational biologist helping a researcher interpret bulk ",
  "RNA-seq differential expression results. Write a clear, rigorous ",
  "biological narrative grounded ONLY in the data provided. Rules: ",
  "(1) Do not invent gene functions or pathways you are not confident about; ",
  "when uncertain, say so. ",
  "(2) Treat the supplied statistics as given; do not claim significance ",
  "beyond what the stated thresholds support. ",
  "(3) Distinguish well-established biology from speculation, and flag ",
  "speculation explicitly. ",
  "(4) Note relevant caveats (correlation is not causation; enrichment ",
  "reflects the ranked/selected gene list, not proof of mechanism; possible ",
  "confounders). ",
  "(5) Be concise and well structured. Output GitHub-flavored Markdown."
)

#' Compact text summary of a DE contrast for the AI prompt
#'
#' @param de a validated DE results data.frame (needs `gene`,
#'   `log2FoldChange`, `padj`; uses `stat` if present)
#' @param top_n max genes listed per direction
#' @param padj_thr,lfc_thr significance thresholds for the up/down counts
#' @return a single character string
#' @keywords internal
summarize_de_for_ai <- function(de, top_n = 30, padj_thr = 0.05, lfc_thr = 1) {
  de <- validate_de_results(de)
  de$gene <- as.character(de$gene)
  ok <- !is.na(de$padj) & !is.na(de$log2FoldChange) & nzchar(de$gene)
  de <- de[ok, , drop = FALSE]
  n_total <- nrow(de)
  sig <- de[de$padj < padj_thr & abs(de$log2FoldChange) > lfc_thr, , drop = FALSE]
  up  <- sig[sig$log2FoldChange > 0, , drop = FALSE]
  dn  <- sig[sig$log2FoldChange < 0, , drop = FALSE]
  up  <- up[order(-up$log2FoldChange), , drop = FALSE]
  dn  <- dn[order(dn$log2FoldChange), , drop = FALSE]
  fmt <- function(d) {
    if (nrow(d) == 0) return("  (none at these thresholds)")
    d <- utils::head(d, top_n)
    paste(sprintf("  %s (log2FC=%+.2f, padj=%.1e)",
                  d$gene, d$log2FoldChange, d$padj), collapse = "\n")
  }
  paste0(
    sprintf(paste0("Genes tested: %d. Significant at padj < %.3g and ",
                   "|log2FC| > %.3g: %d up, %d down.\n\n"),
            n_total, padj_thr, lfc_thr, nrow(up), nrow(dn)),
    sprintf("Top up-regulated genes (max %d):\n%s\n\n", top_n, fmt(up)),
    sprintf("Top down-regulated genes (max %d):\n%s", top_n, fmt(dn))
  )
}

#' Compact text summary of enrichment results for the AI prompt
#'
#' @param enrich an enrichment result list (`method`, `table`) as produced by
#'   the Enrichment tab, or NULL
#' @param n_terms max terms to list
#' @return a character string, or NULL when no usable results are supplied
#' @keywords internal
summarize_enrich_for_ai <- function(enrich, n_terms = 15) {
  if (is.null(enrich) || !is.list(enrich)) return(NULL)
  tab <- enrich$table
  if (is.null(tab) || !is.data.frame(tab) || nrow(tab) == 0) return(NULL)
  method <- enrich$method %||% "gsea"
  tab <- utils::head(tab, n_terms)
  if (identical(method, "gsea")) {
    header <- sprintf(paste0("GSEA -- top %d gene sets by FDR ",
                             "(NES > 0 = enriched on the up-regulated side):"),
                      nrow(tab))
    lines <- sprintf("  %s (NES=%+.2f, padj=%.1e)",
                     tab$pathway, as.numeric(tab$NES), as.numeric(tab$padj))
  } else {
    header <- sprintf("ORA -- top %d over-represented terms by FDR:", nrow(tab))
    desc <- tab$Description %||% tab$ID
    cnt  <- if (!is.null(tab$Count)) tab$Count else rep("?", nrow(tab))
    lines <- sprintf("  %s (genes=%s, padj=%.1e)",
                     desc, cnt, as.numeric(tab$padj))
  }
  paste0(header, "\n", paste(lines, collapse = "\n"))
}

#' Build the interpretation prompt from a contrast (+ optional enrichment)
#'
#' @param de DE results data.frame for the active contrast
#' @param enrich optional enrichment result list (see [summarize_enrich_for_ai()])
#' @param organism organism keyword ("human" / "mouse" / "rat")
#' @param contrast_params optional list with `treated`, `reference`,
#'   `design_var` (used to state the direction of the comparison)
#' @param top_n,n_terms context-size limits passed to the summarizers
#' @return a single character string (the user prompt)
#' @export
build_interpret_prompt <- function(de, enrich = NULL, organism = "human",
                                   contrast_params = NULL, top_n = 30,
                                   n_terms = 15) {
  org <- tryCatch(organism_info(organism)$taxon,
                  error = function(e) organism %||% "the study organism")
  contrast_line <- ""
  if (is.list(contrast_params) && !is.null(contrast_params$treated) &&
      !is.null(contrast_params$reference)) {
    contrast_line <- sprintf(
      paste0("Contrast: %s vs %s (positive log2FC = higher in %s), ",
             "design variable '%s'.\n"),
      contrast_params$treated, contrast_params$reference,
      contrast_params$treated, contrast_params$design_var %||% "group")
  }
  de_block <- summarize_de_for_ai(de, top_n = top_n)
  en_block <- summarize_enrich_for_ai(enrich, n_terms = n_terms)
  en_section <- if (is.null(en_block)) {
    "No functional enrichment results were provided.\n"
  } else {
    paste0("=== Functional enrichment ===\n", en_block, "\n")
  }
  paste0(
    sprintf("Organism: %s.\n", org),
    contrast_line,
    "\n=== Differential expression ===\n", de_block, "\n\n",
    en_section, "\n",
    "Please write a biological interpretation of these results using the ",
    "following Markdown sections:\n",
    "1. **Summary** -- 2-3 sentences on the overall biological signal.\n",
    "2. **Up-regulated programs** -- themes among the up-regulated genes.\n",
    "3. **Down-regulated programs** -- themes among the down-regulated genes.\n",
    "4. **Pathways & enrichment** -- interpret the enriched terms ",
    "(skip this section if no enrichment was provided).\n",
    "5. **Caveats & follow-up** -- statistical/biological caveats and 2-3 ",
    "concrete next experiments.\n",
    "Ground every claim in the genes and terms listed above; flag any ",
    "speculation explicitly."
  )
}

#' Estimate the USD cost of a call from token usage
#'
#' @param input_tokens,output_tokens token counts (may be NA)
#' @param model API model id
#' @return a numeric cost in USD, or NA if it cannot be computed
#' @keywords internal
estimate_cost <- function(input_tokens, output_tokens, model) {
  p <- MODEL_PRICING[[model]]
  if (is.null(p) || is.null(input_tokens) || is.null(output_tokens) ||
      length(input_tokens) != 1 || length(output_tokens) != 1 ||
      is.na(input_tokens) || is.na(output_tokens)) {
    return(NA_real_)
  }
  (input_tokens * p[["in_"]] + output_tokens * p[["out"]]) / 1e6
}

#' Call the Anthropic Claude Messages API
#'
#' A thin \pkg{httr2} wrapper around `POST /v1/messages`. Uses adaptive
#' thinking and returns the assistant's text plus token usage. The API key is
#' never logged or stored; supply it per call.
#'
#' @param prompt the user prompt (character)
#' @param api_key Anthropic API key (defaults to the `ANTHROPIC_API_KEY`
#'   environment variable)
#' @param model API model id (see the AI_MODELS object)
#' @param system optional system prompt
#' @param max_tokens output token cap
#' @param timeout request timeout in seconds
#' @return a list with `text`, `model`, `input_tokens`, `output_tokens`,
#'   `cost_usd`
#' @export
call_claude <- function(prompt, api_key = Sys.getenv("ANTHROPIC_API_KEY"),
                        model = "claude-opus-4-8", system = RNAFLOW_AI_SYSTEM,
                        max_tokens = 6000L, timeout = 120) {
  if (!requireNamespace("httr2", quietly = TRUE)) {
    stop("Package 'httr2' is required for AI interpretation. ",
         "Install with: install.packages('httr2')", call. = FALSE)
  }
  if (!nzchar(api_key %||% "")) {
    stop("No Anthropic API key supplied. Paste your key in the AI tab, or ",
         "set the ANTHROPIC_API_KEY environment variable.", call. = FALSE)
  }

  body <- list(
    model = model,
    max_tokens = as.integer(max_tokens),
    thinking = list(type = "adaptive"),
    messages = list(list(role = "user", content = prompt))
  )
  if (!is.null(system) && nzchar(system)) body$system <- system

  resp <- httr2::request("https://api.anthropic.com/v1/messages") |>
    httr2::req_headers(
      `x-api-key` = api_key,
      `anthropic-version` = "2023-06-01",
      `content-type` = "application/json") |>
    httr2::req_body_json(body) |>
    httr2::req_timeout(timeout) |>
    httr2::req_error(is_error = function(resp) FALSE) |>
    httr2::req_perform()

  status <- httr2::resp_status(resp)
  parsed <- tryCatch(httr2::resp_body_json(resp), error = function(e) NULL)

  if (status != 200) {
    msg <- parsed$error$message %||% httr2::resp_status_desc(resp)
    if (status == 401) {
      msg <- paste0(msg, " (check that your API key is valid)")
    } else if (status == 429) {
      msg <- paste0(msg, " (rate limited -- wait and retry)")
    }
    stop(sprintf("Claude API error (HTTP %d): %s", status, msg), call. = FALSE)
  }
  if (identical(parsed$stop_reason, "refusal")) {
    stop("The model declined to answer this request.", call. = FALSE)
  }

  blocks <- parsed$content %||% list()
  txt <- vapply(blocks, function(b) {
    if (identical(b$type, "text")) b$text %||% "" else ""
  }, character(1))
  text <- paste(txt[nzchar(txt)], collapse = "\n")
  if (!nzchar(text)) {
    stop("The API returned an empty response.", call. = FALSE)
  }

  usage <- parsed$usage %||% list()
  it <- usage$input_tokens  %||% NA_integer_
  ot <- usage$output_tokens %||% NA_integer_
  list(
    text = text,
    model = parsed$model %||% model,
    input_tokens = it,
    output_tokens = ot,
    cost_usd = estimate_cost(it, ot, parsed$model %||% model),
    truncated = identical(parsed$stop_reason, "max_tokens")
  )
}

#' Interpret a DE contrast (+ enrichment) with Claude
#'
#' Convenience wrapper: builds the prompt with [build_interpret_prompt()] and
#' sends it with [call_claude()].
#'
#' @inheritParams build_interpret_prompt
#' @param api_key Anthropic API key
#' @param model API model id
#' @param max_tokens output token cap
#' @param timeout request timeout in seconds
#' @return the [call_claude()] result list, with the `prompt` attached
#' @export
interpret_results <- function(de, enrich = NULL, organism = "human",
                              contrast_params = NULL,
                              api_key = Sys.getenv("ANTHROPIC_API_KEY"),
                              model = "claude-opus-4-8",
                              top_n = 30, n_terms = 15, max_tokens = 6000L,
                              timeout = 120) {
  prompt <- build_interpret_prompt(de, enrich = enrich, organism = organism,
                                   contrast_params = contrast_params,
                                   top_n = top_n, n_terms = n_terms)
  res <- call_claude(prompt, api_key = api_key, model = model,
                     system = RNAFLOW_AI_SYSTEM, max_tokens = max_tokens,
                     timeout = timeout)
  res$prompt <- prompt
  res
}
