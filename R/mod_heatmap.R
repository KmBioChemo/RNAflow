#' Heatmap module
#'
#' @param id namespace ID
#' @name mod_heatmap
NULL

#' @rdname mod_heatmap
#' @export
mod_heatmap_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    ui_page_header(
      "Expression heatmap",
      "Clustered expression of selected genes across samples.",
      about = paste(
        "A heatmap shows whether your genes of interest actually separate the",
        "conditions and cluster into coherent patterns. Row-scaling (z-score",
        "per gene) puts every gene on the same footing, so co-regulated gene",
        "modules and any outlier samples become visible at a glance -- a quick",
        "check that the DE signal is structured, not noise.")),
    sidebar = bslib::sidebar(
      width = 320,
      shiny::radioButtons(ns("src"), "Genes to display",
                          choices = c("Top N by padj" = "top_n",
                                      "All significant" = "all_sig"),
                          selected = "top_n"),
      shiny::conditionalPanel(
        condition = sprintf("input['%s'] == 'top_n'", ns("src")),
        ui_slider_num(ns("n_sld"), ns("n_num"), "Top N genes", 5, 200, 40, 1)
      ),
      ui_slider_num(ns("lfc"), ns("lfc"),
                    "|log2FC| >", 0, 5, 1, 0.1),
      ui_slider_num(ns("padj"), ns("padj"),
                    "padj <", 0.001, 0.2, 0.05, 0.001),
      shiny::selectInput(ns("palette"), "Palette",
                         choices = PALETTE_CHOICES, selected = "RdBu"),
      shiny::checkboxInput(ns("restrict"),
                           "Restrict to active contrast groups", value = FALSE),
      ui_advanced_panel(
        shiny::textInput(ns("title"), "Title", ""),
        shiny::checkboxInput(ns("show_title"), "Show title", TRUE),
        shiny::checkboxInput(ns("cr"), "Cluster rows", TRUE),
        shiny::checkboxInput(ns("cc"), "Cluster columns", TRUE),
        shiny::checkboxInput(ns("rn"), "Show row names", TRUE),
        shiny::checkboxInput(ns("cn"), "Show column names", TRUE),
        shiny::checkboxInput(ns("dr"), "Show row dendrogram", TRUE),
        shiny::checkboxInput(ns("dc"), "Show column dendrogram", TRUE),
        shiny::checkboxInput(ns("lg"), "Show color legend", TRUE),
        shiny::checkboxInput(ns("al"), "Show annotation legend", TRUE),
        shiny::checkboxInput(ns("dir"), "Show direction annotation", FALSE),
        shiny::textInput(ns("ann_title"), "Column annotation header", "")
      ),
      ui_export_bar(ns("hm"), 8, 8)
    ),
    shiny::uiOutput(ns("msg")),
    shinycssloaders::withSpinner(
      shiny::plotOutput(ns("plot"), height = "600px"),
      type = 6, color = "#1D9E75"
    )
  )
}

#' @rdname mod_heatmap
#' @param de_reactive reactive for DE results
#' @param counts_reactive reactive for counts (raw or normalized)
#' @param metadata_reactive reactive for metadata
#' @param contrast_params_reactive optional reactive returning the active
#'   contrast's parameter list (to enable "restrict to contrast")
#' @export
mod_heatmap_server <- function(id, de_reactive, counts_reactive, metadata_reactive,
                               contrast_params_reactive = NULL) {
  shiny::moduleServer(id, function(input, output, session) {

    hm_obj <- shiny::reactive({
      res <- de_reactive()
      shiny::req(res, input$lfc, input$padj)
      counts <- counts_reactive()
      if (is.null(counts)) {
        return(NULL)
      }
      meta <- metadata_reactive()
      if (isTRUE(input$restrict) && !is.null(contrast_params_reactive)) {
        sub <- restrict_to_contrast(counts, meta, contrast_params_reactive())
        counts <- sub$counts; meta <- sub$metadata
      }
      tryCatch({
        fig_heatmap(
          counts_mat = counts, res = res, metadata = meta,
          group_colors = list(),
          n_genes = input$n_num %||% 40,
          gene_src = input$src %||% "top_n",
          padj_thr = input$padj, lfc_thr = input$lfc,
          palette_name = input$palette,
          title = input$title, show_title = input$show_title,
          cluster_rows = input$cr, cluster_cols = input$cc,
          show_rownames = input$rn, show_colnames = input$cn,
          show_dend_rows = input$dr, show_dend_cols = input$dc,
          show_legend = input$lg, annotation_legend = input$al,
          direction_annotation = input$dir,
          ann_title = input$ann_title %||% ""
        )
      }, error = function(e) {
        shiny::showNotification(paste("Heatmap:", conditionMessage(e)),
                                type = "warning", duration = 8)
        NULL
      })
    })

    output$msg <- shiny::renderUI({
      if (is.null(counts_reactive())) {
        ui_banner("Heatmap requires a counts matrix. ",
                  "Upload counts in the Data tab to enable.",
                  type = "warning")
      } else NULL
    })

    output$plot <- shiny::renderPlot({
      ph <- hm_obj()
      shiny::req(ph)
      print(ph)
    }, height = 600, res = 112)

    output$hm_dl <- shiny::downloadHandler(
      filename = function() paste0("heatmap.", input$hm_fmt),
      content  = function(file) {
        ph <- hm_obj(); shiny::req(ph)
        save_pheatmap(ph, file, input$hm_fmt,
                      input$hm_w, input$hm_h, as.integer(input$hm_dpi))
      }
    )
  })
}
