#' Linked explorer module
#'
#' Shiny module rendering a \pkg{crosstalk}-linked volcano and DE table for the
#' active contrast: brushing points in the \pkg{plotly} volcano highlights the
#' matching rows in the table, and the current selection is echoed as a
#' gene list you can copy or download. Wraps the pure [fig_linked] layer.
#'
#' @param id namespace ID
#' @name mod_linked
NULL

#' @rdname mod_linked
#' @export
mod_linked_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 300,
      ui_section_title("Thresholds"),
      ui_slider_num(ns("lfc_sld"), ns("lfc_num"), "|log2FC| >", 0, 5, 1, 0.1),
      ui_slider_num(ns("padj_sld"), ns("padj_num"), "padj <",
                    0.001, 0.2, 0.05, 0.001),
      shiny::tags$hr(style = "margin:8px 0;"),
      shiny::tags$small(
        class = "text-muted",
        "Drag a box (or lasso) on the volcano to select genes -- the table ",
        "filters to your selection, and the list appears below. Click the ",
        "legend to toggle Up / Down / NS."),
      shiny::tags$hr(style = "margin:8px 0;"),
      shiny::downloadButton(ns("dl"), "Download selected genes",
                            class = "btn btn-outline-secondary btn-sm",
                            style = "width:100%;")
    ),
    shiny::uiOutput(ns("msg")),
    shinycssloaders::withSpinner(
      plotly::plotlyOutput(ns("volcano"), height = "420px"),
      type = 6, color = "#1D9E75"
    ),
    shiny::uiOutput(ns("selinfo")),
    shiny::tags$hr(),
    DT::DTOutput(ns("table"))
  )
}

#' @rdname mod_linked
#' @param de_reactive reactive returning the active contrast DE data.frame
#' @export
mod_linked_server <- function(id, de_reactive) {
  shiny::moduleServer(id, function(input, output, session) {

    shiny::observeEvent(input$lfc_sld,
      shiny::updateNumericInput(session, "lfc_num", value = input$lfc_sld),
      ignoreInit = TRUE)
    shiny::observeEvent(input$lfc_num,
      shiny::updateSliderInput(session, "lfc_sld", value = input$lfc_num),
      ignoreInit = TRUE)
    shiny::observeEvent(input$padj_sld,
      shiny::updateNumericInput(session, "padj_num", value = input$padj_sld),
      ignoreInit = TRUE)
    shiny::observeEvent(input$padj_num,
      shiny::updateSliderInput(session, "padj_sld", value = input$padj_num),
      ignoreInit = TRUE)

    # One SharedData drives both widgets; a stable group keeps them linked.
    shared <- shiny::reactive({
      de <- de_reactive()
      shiny::req(de)
      df <- linked_volcano_df(de, padj_thr = input$padj_num %||% 0.05,
                              lfc_thr = input$lfc_num %||% 1)
      crosstalk::SharedData$new(df, key = ~gene, group = session$ns("link"))
    })

    output$msg <- shiny::renderUI({
      if (is.null(de_reactive())) {
        shiny::div(class = "demo-banner",
                   style = "background:#FDEBD0;border-color:#F5B041;color:#7E5109;",
                   "⚠ No active contrast. Run DESeq2 first.")
      } else NULL
    })

    output$volcano <- plotly::renderPlotly({
      fig_linked_volcano(shared())
    })

    # Client-side DT is required for crosstalk linking to work.
    output$table <- DT::renderDT({
      sd <- shared()
      DT::datatable(
        sd, rownames = FALSE, class = "compact stripe hover",
        extensions = "Scroller",
        options = list(pageLength = 12, scrollX = TRUE, deferRender = TRUE))
    }, server = FALSE)

    selected_genes <- shiny::reactive({
      sd <- shared()
      sel <- sd$selection()
      if (is.null(sel) || !any(sel)) return(character(0))
      keys <- sd$key()
      keys[sel]
    })

    output$selinfo <- shiny::renderUI({
      g <- selected_genes()
      if (length(g) == 0) {
        return(shiny::tags$p(
          class = "text-muted", style = "font-size:12px;margin:6px 0;",
          "No genes selected -- brush the volcano to pick some."))
      }
      shiny::tags$p(
        class = "text-muted", style = "font-size:12px;margin:6px 0;",
        sprintf("%d gene%s selected: ", length(g),
                if (length(g) == 1) "" else "s"),
        shiny::tags$code(paste(utils::head(g, 40), collapse = ", ")),
        if (length(g) > 40) sprintf(" ... (+%d)", length(g) - 40) else NULL)
    })

    output$dl <- shiny::downloadHandler(
      filename = function() "selected_genes.txt",
      content = function(file) {
        g <- selected_genes()
        writeLines(if (length(g)) g else "# no genes selected", file)
      }
    )
  })
}
