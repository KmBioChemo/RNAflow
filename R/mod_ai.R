#' AI interpretation module
#'
#' Shiny module wrapping the pure [ai_interpret] layer. Sends a compact summary
#' of the active contrast (top gene names + fold-changes, and the latest
#' enrichment terms) to Anthropic's Claude API and renders the returned
#' biological narrative. The API key lives only in the session -- it is read
#' from a password field or the `ANTHROPIC_API_KEY` environment variable, and
#' is never written to disk.
#'
#' @param id namespace ID
#' @name mod_ai
NULL

#' @rdname mod_ai
#' @export
mod_ai_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 340,
      ui_section_title("Claude API"),
      shiny::passwordInput(ns("api_key"), "API key (session only)",
                           placeholder = "sk-ant-..."),
      shiny::tags$small(
        class = "text-muted",
        "Pasted here it stays in this session only -- never saved to disk. ",
        "Or set the ", shiny::tags$code("ANTHROPIC_API_KEY"),
        " environment variable."),
      shiny::selectInput(ns("model"), "Model", choices = AI_MODELS,
                         selected = "claude-opus-4-8"),
      shiny::tags$hr(style = "margin:8px 0;"),

      ui_section_title("Context sent"),
      ui_slider_num(ns("topn_sld"), ns("topn_num"),
                    "Top genes / direction", 10, 60, 30, 5),
      ui_slider_num(ns("nterm_sld"), ns("nterm_num"),
                    "Top enrichment terms", 5, 30, 15, 1),
      shiny::checkboxInput(ns("use_enrich"),
                           "Include enrichment results (if run)", TRUE),
      shiny::tags$hr(style = "margin:8px 0;"),

      shiny::actionButton(ns("run"), "Interpret with AI",
                          class = "btn btn-primary btn-sm",
                          style = "width:100%;",
                          icon = shiny::icon("wand-magic-sparkles")),
      shiny::tags$hr(style = "margin:8px 0;"),
      shiny::downloadButton(ns("dl"), "Download .md",
                            class = "btn btn-outline-secondary btn-sm",
                            style = "width:100%;")
    ),
    shiny::div(
      class = "demo-banner", style = "margin-bottom:8px;",
      "\u2139 This sends a compact summary of your active contrast (top gene ",
      "names, fold-changes, and enrichment terms) to Anthropic's Claude API. ",
      "Your count matrix and sample metadata never leave your machine. ",
      "AI-generated text can be wrong -- treat it as a hypothesis-generating ",
      "draft, not a conclusion."),
    shiny::uiOutput(ns("notice")),
    shinycssloaders::withSpinner(
      shiny::uiOutput(ns("out")),
      type = 6, color = "#1D9E75"
    ),
    shiny::uiOutput(ns("meta"))
  )
}

#' @rdname mod_ai
#' @param de_reactive reactive returning the active contrast DE data.frame
#' @param enrich_reactive optional reactive returning the latest enrichment
#'   result list (`method`, `table`), as exposed by [mod_enrich_server()]
#' @param organism_reactive optional reactive returning the organism keyword
#' @param contrast_params_reactive optional reactive returning the active
#'   contrast's parameter list (design_var / treated / reference)
#' @param settings_store optional `reactiveVal` holding a settings list; the
#'   latest interpretation is recorded under `$ai_interpretation` so it can be
#'   archived in the HTML report
#' @export
mod_ai_server <- function(id, de_reactive, enrich_reactive = NULL,
                          organism_reactive = NULL,
                          contrast_params_reactive = NULL,
                          settings_store = NULL) {
  shiny::moduleServer(id, function(input, output, session) {

    mirror <- function(s, n) {
      shiny::observeEvent(input[[s]],
        shiny::updateNumericInput(session, n, value = input[[s]]),
        ignoreInit = TRUE)
      shiny::observeEvent(input[[n]],
        shiny::updateSliderInput(session, s, value = input[[n]]),
        ignoreInit = TRUE)
    }
    mirror("topn_sld", "topn_num")
    mirror("nterm_sld", "nterm_num")

    result <- shiny::reactiveVal(NULL)  # list from interpret_results()

    resolve_key <- function() {
      k <- input$api_key %||% ""
      if (!nzchar(k)) k <- Sys.getenv("ANTHROPIC_API_KEY")
      k
    }

    shiny::observeEvent(input$run, {
      de <- de_reactive()
      if (is.null(de)) {
        shiny::showNotification("No active contrast. Run DESeq2 first.",
                                type = "warning"); return()
      }
      key <- resolve_key()
      if (!nzchar(key)) {
        shiny::showNotification(
          "Paste your Anthropic API key (or set ANTHROPIC_API_KEY).",
          type = "warning"); return()
      }
      org <- if (!is.null(organism_reactive)) {
        organism_reactive() %||% "human"
      } else "human"
      # Match the Enrichment tab: convert Ensembl / ENTREZ IDs to gene symbols
      # so the model sees readable gene names (no-op if already symbols).
      de <- tryCatch(map_de_to_symbols(de, org), error = function(e) de)
      enr <- if (isTRUE(input$use_enrich) && !is.null(enrich_reactive)) {
        tryCatch(enrich_reactive(), error = function(e) NULL)
      } else NULL
      cp <- if (!is.null(contrast_params_reactive)) {
        tryCatch(contrast_params_reactive(), error = function(e) NULL)
      } else NULL

      shiny::withProgress(message = "Asking Claude...", value = 0.4, {
        tryCatch({
          res <- interpret_results(
            de, enrich = enr, organism = org, contrast_params = cp,
            api_key = key, model = input$model %||% "claude-opus-4-8",
            top_n = as.integer(input$topn_num %||% 30),
            n_terms = as.integer(input$nterm_num %||% 15))
          result(res)
          if (!is.null(settings_store)) {
            s <- settings_store()
            # Store provenance alongside the text so a saved project / report
            # records how the interpretation was produced. The API key is
            # never included.
            s$ai_interpretation <- list(
              text = res$text, model = res$model,
              generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
              top_n = as.integer(input$topn_num %||% 30),
              n_terms = as.integer(input$nterm_num %||% 15),
              use_enrich = isTRUE(input$use_enrich) && !is.null(enr),
              input_tokens = res$input_tokens,
              output_tokens = res$output_tokens,
              cost_usd = res$cost_usd)
            settings_store(s)
          }
        }, error = function(e) {
          result(NULL)
          shiny::showNotification(
            paste("AI interpretation failed:", conditionMessage(e)),
            type = "error", duration = 10)
        })
      })
    })

    output$notice <- shiny::renderUI({
      if (!is.null(result())) return(NULL)
      shiny::div(class = "demo-banner", style = "margin-bottom:8px;",
                 "Paste your API key, then click ",
                 shiny::strong("Interpret with AI"),
                 " to summarize the active contrast.")
    })

    output$out <- shiny::renderUI({
      res <- result()
      shiny::req(res)
      shiny::div(
        class = "card",
        style = "padding:16px 20px;border:1px solid #e9ecef;border-radius:8px;",
        shiny::markdown(res$text)
      )
    })

    output$meta <- shiny::renderUI({
      res <- result()
      shiny::req(res)
      cost <- if (is.na(res$cost_usd)) "n/a" else sprintf("$%.4f", res$cost_usd)
      trunc_note <- if (isTRUE(res$truncated)) {
        shiny::span(style = "color:#C0392B;",
                    " \u26A0 response hit the token limit and may be cut off.")
      } else NULL
      shiny::tags$p(
        class = "text-muted", style = "font-size:12px;margin-top:8px;",
        sprintf("Model: %s | tokens: %s in / %s out | est. cost: %s",
                res$model,
                res$input_tokens %||% "?", res$output_tokens %||% "?", cost),
        trunc_note)
    })

    output$dl <- shiny::downloadHandler(
      filename = function() "rnaflow_ai_interpretation.md",
      content = function(file) {
        res <- result()
        if (is.null(res)) {
          writeLines("No interpretation generated yet.", file); return()
        }
        writeLines(c(
          "# RNAflow -- AI interpretation",
          sprintf("_Generated with %s via the Anthropic Claude API._", res$model),
          "", res$text), file)
      }
    )

    invisible(result)
  })
}
