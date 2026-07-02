#' Multi-contrast comparison module
#'
#' Shiny module that consumes the contrast store and exposes the four
#' cross-contrast views: Venn, UpSet, side-by-side volcano grid, and the
#' log2FoldChange signature heatmap. Wraps the pure [fig_compare] functions.
#'
#' @param id namespace ID
#' @name mod_compare
NULL

#' @rdname mod_compare
#' @export
mod_compare_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    ui_page_header("Compare contrasts",
                   "Overlap and log2 fold-change across saved DE contrasts."),
    sidebar = bslib::sidebar(
      width = 330,
      ui_section_title("Contrasts to compare"),
      shiny::uiOutput(ns("contrast_picker")),
      shiny::tags$hr(style = "margin:8px 0;"),
      ui_section_title("View"),
      shiny::radioButtons(
        ns("view"), NULL,
        choices = c("Venn diagram"      = "venn",
                    "UpSet plot"        = "upset",
                    "Volcano grid"      = "grid",
                    "log2FC heatmap"    = "heatmap"),
        selected = "venn"
      ),
      shiny::tags$hr(style = "margin:8px 0;"),
      ui_section_title("Significance"),
      ui_slider_num(ns("lfc_sld"), ns("lfc_num"),
                    "|log2 fold change| >", 0, 5, 1, 0.1),
      ui_slider_num(ns("padj_sld"), ns("padj_num"),
                    "FDR (padj) <", 0.001, 0.2, 0.05, 0.001),
      shiny::radioButtons(
        ns("direction"), "Direction",
        choices = c("Either" = "either", "Up only" = "up", "Down only" = "down"),
        selected = "either", inline = TRUE
      ),
      # View-specific controls
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'upset'", ns("view")),
        shiny::numericInput(ns("min_size"), "Min intersection size",
                            value = 1, min = 1, step = 1)
      ),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'grid'", ns("view")),
        ui_slider_num(ns("nlab_sld"), ns("nlab_num"),
                      "Top genes to label", 0, 30, 8, 1),
        shiny::checkboxInput(ns("publication"),
                             "Publication mode (export)", value = FALSE)
      ),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'heatmap'", ns("view")),
        shiny::radioButtons(
          ns("gene_src"), "Genes",
          choices = c("Significant union" = "sig_union",
                      "Most variable"     = "top_var"),
          selected = "sig_union"
        ),
        ui_slider_num(ns("ngenes_sld"), ns("ngenes_num"),
                      "Max genes", 5, 200, 50, 5),
        shiny::selectInput(ns("palette"), "Palette",
                           choices = c("RdBu", "RdYlBu", "PuOr", "BrBG",
                                       "PRGn", "Spectral"),
                           selected = "RdBu")
      ),
      ui_export_bar(ns("cmp"), 9, 6)
    ),
    shiny::uiOutput(ns("notice")),
    shinycssloaders::withSpinner(
      shiny::plotOutput(ns("plot"), height = "600px"),
      type = 6, color = "#1D9E75"
    )
  )
}

#' @rdname mod_compare
#' @param store_reactive a reactive returning the contrast store (named list,
#'   as built by [contrast_store_upsert()])
#' @export
mod_compare_server <- function(id, store_reactive) {
  shiny::moduleServer(id, function(input, output, session) {

    mirror <- function(s, n) {
      shiny::observeEvent(input[[s]],
        shiny::updateNumericInput(session, n, value = input[[s]]),
        ignoreInit = TRUE)
      shiny::observeEvent(input[[n]],
        shiny::updateSliderInput(session, s, value = input[[n]]),
        ignoreInit = TRUE)
    }
    mirror("lfc_sld", "lfc_num")
    mirror("padj_sld", "padj_num")
    mirror("nlab_sld", "nlab_num")
    mirror("ngenes_sld", "ngenes_num")

    # Contrast multi-select, rebuilt as the store grows
    output$contrast_picker <- shiny::renderUI({
      store <- store_reactive()
      ns <- session$ns
      if (length(store) == 0) {
        return(shiny::div(class = "demo-banner",
                          "No contrasts yet. Run DESeq2 on the Data tab; ",
                          "each run is saved here automatically."))
      }
      labels <- names(store)
      shiny::checkboxGroupInput(
        ns("contrasts"), NULL,
        choices  = labels,
        selected = utils::head(labels, 2)
      )
    })

    # Named list of DE data.frames for the current selection
    selected_contrasts <- shiny::reactive({
      store <- store_reactive()
      sel <- input$contrasts
      sel <- sel[sel %in% names(store)]
      contrast_store_results(store[sel])
    })

    # Build the requested figure, surfacing problems as clean notices
    build <- shiny::reactive({
      contrasts <- selected_contrasts()
      view <- input$view %||% "venn"
      min_needed <- 2L
      if (length(contrasts) < min_needed) {
        return(list(error = sprintf(
          "Select at least %d contrasts to compare.", min_needed)))
      }
      shiny::req(input$lfc_num, input$padj_num)
      obj <- tryCatch({
        switch(
          view,
          venn = fig_venn(
            contrast_sig_sets(contrasts, input$padj_num, input$lfc_num,
                              input$direction)),
          upset = fig_upset(
            contrast_sig_sets(contrasts, input$padj_num, input$lfc_num,
                              input$direction),
            min_size = max(1L, as.integer(input$min_size %||% 1))),
          grid = fig_volcano_grid(
            contrasts, lfc_thr = input$lfc_num, padj_thr = input$padj_num,
            n_label = max(0L, as.integer(input$nlab_num %||% 8)),
            mode = if (isTRUE(input$publication)) "publication" else "exploration"),
          heatmap = fig_lfc_heatmap(
            contrasts, gene_src = input$gene_src %||% "sig_union",
            n_genes = max(2L, as.integer(input$ngenes_num %||% 50)),
            padj_thr = input$padj_num, lfc_thr = input$lfc_num,
            palette_name = input$palette %||% "RdBu")
        )
      }, error = function(e) structure(conditionMessage(e), class = "cmp_error"))
      if (inherits(obj, "cmp_error")) return(list(error = as.character(obj)))
      list(obj = obj)
    })

    output$notice <- shiny::renderUI({
      res <- build()
      if (is.null(res$error)) return(NULL)
      shiny::div(class = "demo-banner",
                 style = "margin-bottom:8px;", res$error)
    })

    output$plot <- shiny::renderPlot({
      res <- build()
      shiny::validate(shiny::need(is.null(res$error), res$error))
      draw_compare(res$obj)
    })

    output$cmp_dl <- shiny::downloadHandler(
      filename = function() paste0("contrast_comparison.", input$cmp_fmt),
      content  = function(file) {
        res <- build()
        if (!is.null(res$error)) {
          stop("Cannot export: ", res$error, call. = FALSE)
        }
        save_compare(res$obj, file, input$cmp_fmt,
                     input$cmp_w, input$cmp_h, as.integer(input$cmp_dpi))
      }
    )
  })
}
