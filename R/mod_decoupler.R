#' Activity inference module
#'
#' Shiny module wrapping the pure [analysis_decoupler] layer. Infers
#' transcription-factor (CollecTRI) or pathway (PROGENy) activity for the
#' active contrast with \pkg{decoupleR} and renders a diverging bar chart plus
#' a results table.
#'
#' @param id namespace ID
#' @name mod_activity
NULL

#' @rdname mod_activity
#' @export
mod_activity_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 340,
      shiny::uiOutput(ns("organism_note")),
      shiny::radioButtons(
        ns("type"), "Activity",
        choices = c("Transcription factors (CollecTRI)" = "tf",
                    "Pathways (PROGENy)" = "pathway"),
        selected = "tf"),
      shiny::selectInput(ns("rank_by"), "Rank genes by",
                         choices = c("Wald statistic" = "stat",
                                     "signed -log10(p)" = "signed_p",
                                     "log2 fold change" = "log2fc"),
                         selected = "stat"),
      shiny::tags$small(
        class = "text-muted",
        "Prior-knowledge networks are fetched from OmniPath ",
        "(needs internet on first use; then cached for the session)."),
      shiny::tags$hr(style = "margin:8px 0;"),
      shiny::actionButton(ns("run"), "Infer activity",
                          class = "btn btn-primary btn-sm",
                          style = "width:100%;"),
      shiny::tags$hr(style = "margin:8px 0;"),
      ui_slider_num(ns("n_sld"), ns("n_num"), "Top terms", 5, 40, 20, 1),
      shiny::checkboxInput(ns("publication"), "Publication mode (export)",
                           FALSE),
      ui_export_bar(ns("act"), 7, 6)
    ),
    shiny::uiOutput(ns("notice")),
    shinycssloaders::withSpinner(
      shiny::plotOutput(ns("plot"), height = "540px"),
      type = 6, color = "#1D9E75"
    ),
    shiny::tags$hr(),
    DT::DTOutput(ns("table"))
  )
}

#' @rdname mod_activity
#' @param de_reactive reactive returning the active contrast DE data.frame
#' @param organism_reactive reactive returning the organism keyword
#' @export
mod_activity_server <- function(id, de_reactive, organism_reactive) {
  shiny::moduleServer(id, function(input, output, session) {

    shiny::observeEvent(input$n_sld,
      shiny::updateNumericInput(session, "n_num", value = input$n_sld),
      ignoreInit = TRUE)
    shiny::observeEvent(input$n_num,
      shiny::updateSliderInput(session, "n_sld", value = input$n_num),
      ignoreInit = TRUE)

    # Cache fetched networks per (type, organism) for the session.
    net_cache <- shiny::reactiveVal(list())
    result <- shiny::reactiveVal(NULL)   # list(type, table)

    output$organism_note <- shiny::renderUI({
      org <- organism_reactive() %||% "human"
      shiny::div(class = "demo-banner", style = "margin-bottom:6px;",
                 sprintf("Organism: %s (set on the Data tab)", org))
    })

    fetch_network <- function(type, org) {
      key <- paste(type, org, sep = ":")
      cache <- net_cache()
      if (!is.null(cache[[key]])) return(cache[[key]])
      net <- if (type == "tf") get_tf_network(org) else get_pathway_network(org)
      cache[[key]] <- net
      net_cache(cache)
      net
    }

    shiny::observeEvent(input$run, {
      de <- de_reactive()
      org <- organism_reactive() %||% "human"
      if (is.null(de)) {
        shiny::showNotification("No active contrast. Run DESeq2 first.",
                                type = "warning"); return()
      }
      # Activity inference works on gene symbols, like enrichment.
      de <- tryCatch(map_de_to_symbols(de, org), error = function(e) de)
      shiny::withProgress(message = "Fetching network + inferring activity...",
                          value = 0.3, {
        tryCatch({
          net <- fetch_network(input$type, org)
          shiny::incProgress(0.4)
          meth <- if (input$type == "tf") "ulm" else "mlm"
          mcol <- if (input$type == "tf") "mor" else "weight"
          tab <- run_activity(de, net, method = meth, mor_col = mcol,
                              by = input$rank_by)
          result(list(type = input$type, table = tab))
          shiny::showNotification(
            sprintf("Activity inferred: %d %s.", nrow(tab),
                    if (input$type == "tf") "regulators" else "pathways"),
            type = "message", duration = 4)
        }, error = function(e) {
          result(NULL)
          shiny::showNotification(
            paste("Activity inference failed:", conditionMessage(e)),
            type = "error", duration = 10)
        })
      })
    })

    output$notice <- shiny::renderUI({
      r <- result()
      if (is.null(r)) {
        return(shiny::div(class = "demo-banner", style = "margin-bottom:8px;",
                          "Pick an activity type and click Infer activity. ",
                          "TF activity uses CollecTRI regulons; pathway ",
                          "activity uses PROGENy."))
      }
      if (nrow(r$table) == 0) {
        return(shiny::div(class = "demo-banner", style = "margin-bottom:8px;",
                          "No activity scores returned."))
      }
      NULL
    })

    cur_plot <- shiny::reactive({
      r <- result()
      shiny::req(r, nrow(r$table) > 0)
      mode <- if (isTRUE(input$publication)) "publication" else "exploration"
      ttl <- if (r$type == "tf") "Transcription-factor activity" else
        "Pathway activity"
      fig_activity_bar(r$table, n = max(5L, as.integer(input$n_num %||% 20)),
                       title = ttl, mode = mode)
    })

    output$plot <- shiny::renderPlot({ print(cur_plot()) })

    output$table <- DT::renderDT({
      r <- result()
      shiny::req(r, nrow(r$table) > 0)
      df <- r$table
      df$score <- round(df$score, 3)
      for (col in c("p_value", "padj")) {
        if (col %in% colnames(df)) df[[col]] <- formatC(df[[col]], format = "e",
                                                        digits = 2)
      }
      DT::datatable(df, rownames = FALSE, class = "compact stripe hover",
                    options = list(pageLength = 10, scrollX = TRUE))
    })

    output$act_dl <- shiny::downloadHandler(
      filename = function() paste0("activity_", result()$type %||% "plot",
                                   ".", input$act_fmt),
      content = function(file) {
        save_ggplot(cur_plot(), file, input$act_fmt,
                    input$act_w, input$act_h, as.integer(input$act_dpi))
      }
    )
  })
}
