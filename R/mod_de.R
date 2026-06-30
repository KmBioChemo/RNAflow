#' Differential expression module
#'
#' UI + server for running DESeq2 from counts + metadata. Lets the user
#' pick the design variable, contrast levels, and shrinkage options.
#' Returns a reactive holding the tidy results data.frame.
#'
#' @param id namespace ID
#' @name mod_de
NULL

#' @rdname mod_de
#' @export
mod_de_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::card(
    bslib::card_header("Differential expression (DESeq2)"),
    bslib::card_body(
      shiny::uiOutput(ns("design_ui")),
      shiny::uiOutput(ns("contrast_ui")),
      ui_advanced_panel(
        shiny::checkboxInput(ns("shrink"), "Apply LFC shrinkage (apeglm)",
                             value = TRUE),
        shiny::numericInput(ns("min_count"), "Minimum row sum to keep gene",
                            value = 10, min = 0, step = 1),
        shiny::numericInput(ns("alpha"), "FDR threshold (independent filtering)",
                            value = 0.05, min = 0.001, max = 0.5, step = 0.01)
      ),
      shiny::actionButton(ns("run"), "Run DESeq2",
                          class = "btn btn-primary btn-sm",
                          style = "margin-top:8px;"),
      shiny::uiOutput(ns("status"))
    )
  )
}

#' @rdname mod_de
#' @param data_mod the value returned by [mod_data_server()]
#' @export
mod_de_server <- function(id, data_mod) {
  shiny::moduleServer(id, function(input, output, session) {

    de_r <- shiny::reactiveVal(NULL)

    # Build design dropdown from metadata columns
    output$design_ui <- shiny::renderUI({
      meta <- data_mod$metadata()
      if (is.null(meta) || ncol(meta) < 2) {
        return(shiny::div(class = "demo-banner",
                          "Load counts + metadata to enable DE."))
      }
      ns <- session$ns
      choices <- colnames(meta)[-1]
      shiny::selectInput(ns("design_var"),
                         "Design variable (last term of ~ X)",
                         choices = choices, selected = choices[1])
    })

    output$contrast_ui <- shiny::renderUI({
      meta <- data_mod$metadata()
      shiny::req(meta, input$design_var)
      lv <- unique(as.character(meta[[input$design_var]]))
      if (length(lv) < 2) {
        return(shiny::div(class = "demo-banner",
                          "Variable has fewer than 2 levels."))
      }
      ns <- session$ns
      shiny::tagList(
        shiny::selectInput(ns("level_treated"), "Treated (numerator)",
                           choices = lv, selected = lv[2]),
        shiny::selectInput(ns("level_reference"), "Reference (denominator)",
                           choices = lv, selected = lv[1])
      )
    })

    shiny::observeEvent(input$run, {
      counts <- data_mod$counts()
      meta   <- data_mod$metadata()
      shiny::req(counts, meta, input$design_var,
                 input$level_treated, input$level_reference)
      if (input$level_treated == input$level_reference) {
        shiny::showNotification("Treated and reference levels are the same.",
                                type = "warning")
        return()
      }
      shiny::withProgress(message = "Running DESeq2...", value = 0.3, {
        tryCatch({
          res <- run_deseq2(
            counts, meta,
            design   = stats::as.formula(paste0("~", input$design_var)),
            contrast = c(input$design_var, input$level_treated, input$level_reference),
            shrink   = isTRUE(input$shrink),
            min_count = as.integer(input$min_count %||% 10),
            alpha    = as.numeric(input$alpha %||% 0.05)
          )
          shiny::incProgress(0.7)
          de_r(res)
          shiny::showNotification(
            sprintf("DESeq2 done: %d genes (%d sig at padj < %.2f)",
                    nrow(res),
                    sum(res$padj < (input$alpha %||% 0.05), na.rm = TRUE),
                    input$alpha %||% 0.05),
            type = "message", duration = 5
          )
        }, error = function(e) {
          shiny::showNotification(paste("DESeq2 failed:", conditionMessage(e)),
                                  type = "error", duration = 12)
        })
      })
    })

    output$status <- shiny::renderUI({
      d <- de_r()
      if (is.null(d)) return(NULL)
      shiny::div(class = "demo-banner",
                 sprintf("\u2713 %d genes in results table", nrow(d)))
    })

    list(de_results = shiny::reactive(de_r()))
  })
}
