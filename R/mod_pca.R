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
    ui_page_header(
      "Sample overview (PCA / UMAP)",
      "Low-dimensional map of samples from the whole expression matrix.",
      about = paste(
        "PCA and UMAP compress thousands of genes into a 2-3D map of your",
        "samples, exposing the dominant axes of variation. If replicates",
        "cluster tightly and conditions separate, the design is sound; a",
        "mislabelled sample or an unwanted batch effect usually shows up here",
        "first, before it contaminates the differential-expression results.")),
    sidebar = bslib::sidebar(
      width = 280,
      shiny::selectInput(ns("method"), "Embedding",
                         c("PCA (2D)" = "pca", "PCA (3D)" = "pca3d",
                           "UMAP" = "umap")),
      ui_slider_num(ns("n_sld"), ns("n_num"),
                    "Top variable genes", 50, 5000, 500, 50),
      shiny::uiOutput(ns("col_ui")),
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'umap'", ns("method")),
        ui_slider_num(ns("nn_sld"), ns("nn_num"),
                      "UMAP neighbors", 2, 50, 15, 1),
        shiny::sliderInput(ns("min_dist"), "UMAP min. distance",
                           min = 0.01, max = 0.99, value = 0.1, step = 0.01)
      ),
      shiny::checkboxInput(ns("labels"), "Show sample labels", value = TRUE),
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
        "Load a counts matrix to generate the sample overview."
      ))
      shiny::req(input$n_num)
      method <- input$method %||% "pca"
      if (method %in% c("umap", "pca3d")) {
        shiny::validate(shiny::need(
          ncol(d$counts) >= 4,
          "UMAP and 3D PCA need at least 4 samples."))
      }
      col_by <- if (!is.null(input$col_by) && input$col_by != "(none)") input$col_by else NULL
      labs <- isTRUE(input$labels)
      switch(
        method,
        pca3d = fig_pca_3d(d$counts, d$metadata, input$n_num, col_by,
                           input$title, show_labels = labs),
        umap  = fig_umap(d$counts, d$metadata, input$n_num,
                         n_neighbors = input$nn_num %||% 15,
                         min_dist = input$min_dist %||% 0.1,
                         color_by = col_by, title = input$title,
                         show_labels = labs),
        fig_pca(d$counts, d$metadata, input$n_num, col_by, input$title,
                show_labels = labs)
      )
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
