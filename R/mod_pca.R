#' PCA module
#'
#' @param id namespace ID
#' @name mod_pca
NULL

#' @rdname mod_pca
#' @export
mod_pca_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 280,
      ui_slider_num(ns("n_sld"), ns("n_num"),
                    "Top variable genes", 50, 5000, 500, 50),
      shiny::uiOutput(ns("col_ui")),
      shiny::checkboxInput(ns("restrict"),
                           "Restrict to active contrast groups", value = FALSE),
      shiny::textInput(ns("title"), "Title", ""),
      ui_export_bar(ns("pca"), 7, 6)
    ),
    shiny::uiOutput(ns("msg")),
    shinycssloaders::withSpinner(
      plotly::plotlyOutput(ns("plot"), height = "560px"),
      type = 6, color = "#1D9E75"
    )
  )
}

#' @rdname mod_pca
#' @param counts_reactive reactive for counts matrix (normalized recommended)
#' @param metadata_reactive reactive for metadata
#' @param contrast_params_reactive optional reactive returning the active
#'   contrast's parameter list (to enable "restrict to contrast")
#' @export
mod_pca_server <- function(id, counts_reactive, metadata_reactive,
                           contrast_params_reactive = NULL) {
  shiny::moduleServer(id, function(input, output, session) {

    # counts + metadata, optionally restricted to the active contrast's groups
    data_r <- shiny::reactive({
      counts <- counts_reactive(); meta <- metadata_reactive()
      if (isTRUE(input$restrict) && !is.null(contrast_params_reactive)) {
        sub <- restrict_to_contrast(counts, meta, contrast_params_reactive())
        return(sub)
      }
      list(counts = counts, metadata = meta)
    })

    output$col_ui <- shiny::renderUI({
      meta <- metadata_reactive()
      ns <- session$ns
      ch <- if (!is.null(meta) && ncol(meta) >= 2) {
        c("(none)", colnames(meta)[-1])
      } else "(none)"
      sel <- if (length(ch) > 1) ch[2] else "(none)"
      shiny::selectInput(ns("col_by"), "Color samples by",
                         choices = ch, selected = sel)
    })

    pca_plot <- shiny::reactive({
      d <- data_r()
      shiny::validate(shiny::need(
        !is.null(d$counts),
        "Load a counts matrix to generate the PCA."
      ))
      shiny::req(input$n_num)
      col_by <- if (!is.null(input$col_by) && input$col_by != "(none)") input$col_by else NULL
      fig_pca(d$counts, d$metadata, input$n_num, col_by, input$title)
    })

    output$msg <- shiny::renderUI({
      if (is.null(counts_reactive())) {
        ui_banner("PCA requires a counts matrix.", type = "warning")
      } else NULL
    })

    output$plot <- plotly::renderPlotly({ pca_plot() })

    output$pca_dl <- shiny::downloadHandler(
      filename = function() paste0("pca_plot.", input$pca_fmt),
      content  = function(file) {
        p <- pca_plot()
        fmt <- input$pca_fmt
        if (fmt == "pdf") {
          tf <- tempfile(fileext = ".html")
          htmlwidgets::saveWidget(p, tf, selfcontained = TRUE)
          shiny::showNotification(
            "PDF export: use the camera icon in the PCA toolbar to download PNG.",
            type = "message", duration = 6)
          file.copy(tf, file)
        } else {
          tryCatch(
            plotly::save_image(p, file,
                               width = as.integer(input$pca_w * 96),
                               height = as.integer(input$pca_h * 96), scale = 3),
            error = function(e) {
              htmlwidgets::saveWidget(p, file, selfcontained = TRUE)
              shiny::showNotification(
                "Static export unavailable. Saved as HTML instead.",
                type = "warning", duration = 8)
            }
          )
        }
      }
    )
  })
}
