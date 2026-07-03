#' QC / diagnostics module
#'
#' Shiny module exposing the [fig_qc] diagnostics: p-value histogram, MA plot,
#' sample-correlation heatmap, and library sizes. Helps sanity-check a run
#' before interpreting the results.
#'
#' @param id namespace ID
#' @name mod_qc
NULL

#' @rdname mod_qc
#' @export
mod_qc_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    ui_page_header("Quality control",
                   "Sample correlations, library sizes, and p-value diagnostics."),
    sidebar = bslib::sidebar(
      width = 300,
      ui_section_title("Diagnostic"),
      shiny::radioButtons(
        ns("view"), NULL,
        choices = c("P-value histogram" = "pval",
                    "MA plot"            = "ma",
                    "Sample correlation" = "cor",
                    "Library sizes"      = "lib",
                    "Gene expression"    = "gene"),
        selected = "pval"),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'ma'", ns("view")),
        ui_slider_num(ns("padj_sld"), ns("padj_num"), "FDR (padj) <",
                      0.001, 0.2, 0.05, 0.001)),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'cor'", ns("view")),
        shiny::radioButtons(ns("cor_method"), "Correlation",
                            choices = c("pearson", "spearman"),
                            selected = "pearson", inline = TRUE)),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'gene'", ns("view")),
        shiny::selectizeInput(ns("gene"), "Gene", choices = NULL,
                              options = list(placeholder = "type a gene...")),
        shiny::radioButtons(ns("gene_style"), "Style",
                            choices = c("Raincloud" = "raincloud",
                                        "Beeswarm" = "beeswarm", "Box" = "box"),
                            selected = "raincloud", inline = TRUE)),
      ui_export_bar(ns("qc"), 7, 5)
    ),
    shiny::uiOutput(ns("notice")),
    shinycssloaders::withSpinner(
      shiny::plotOutput(ns("plot"), height = "520px"),
      type = 6, color = "#1D9E75")
  )
}

#' @rdname mod_qc
#' @param de_reactive reactive: active contrast DE data.frame
#' @param counts_reactive reactive: raw counts matrix
#' @param counts_norm_reactive reactive: normalized counts matrix
#' @param metadata_reactive reactive: sample metadata
#' @export
mod_qc_server <- function(id, de_reactive, counts_reactive,
                          counts_norm_reactive, metadata_reactive) {
  shiny::moduleServer(id, function(input, output, session) {

    shiny::observeEvent(input$padj_sld,
      shiny::updateNumericInput(session, "padj_num", value = input$padj_sld),
      ignoreInit = TRUE)
    shiny::observeEvent(input$padj_num,
      shiny::updateSliderInput(session, "padj_sld", value = input$padj_num),
      ignoreInit = TRUE)

    needs <- list(
      pval = "Run DESeq2 to see the p-value distribution.",
      ma   = "Run DESeq2 to see the MA plot.",
      cor  = "Load counts (and run DESeq2 for normalization) for the correlation heatmap.",
      lib  = "Load a counts matrix to see library sizes.",
      gene = "Load counts (normalized) to plot a gene's expression.")

    # Populate the gene selector (server-side for large matrices).
    shiny::observeEvent(counts_norm_reactive(), {
      cn <- counts_norm_reactive()
      if (is.null(cn) || is.null(rownames(cn))) return()
      shiny::updateSelectizeInput(
        session, "gene", choices = rownames(cn), server = TRUE,
        selected = (shiny::isolate(input$gene) %||% rownames(cn)[1]))
    }, ignoreNULL = FALSE)

    cur_plot <- shiny::reactive({
      v <- input$view %||% "pval"
      mode <- "exploration"
      switch(
        v,
        pval = { res <- de_reactive()
                 shiny::validate(shiny::need(!is.null(res), needs$pval))
                 fig_pval_hist(res, mode = mode) },
        ma   = { res <- de_reactive()
                 shiny::validate(shiny::need(!is.null(res), needs$ma))
                 fig_ma(res, padj_thr = input$padj_num %||% 0.05, mode = mode) },
        cor  = { cn <- counts_norm_reactive()
                 shiny::validate(shiny::need(!is.null(cn), needs$cor))
                 fig_sample_cor(cn, metadata_reactive(),
                                method = input$cor_method %||% "pearson") },
        lib  = { ct <- counts_reactive()
                 shiny::validate(shiny::need(!is.null(ct), needs$lib))
                 fig_lib_sizes(ct, metadata_reactive(), mode = mode) },
        gene = { cn <- counts_norm_reactive()
                 shiny::validate(shiny::need(!is.null(cn), needs$gene))
                 shiny::req(input$gene)
                 shiny::validate(shiny::need(
                   input$gene %in% rownames(cn),
                   "Pick a gene present in the counts matrix."))
                 fig_gene_expression(cn, metadata_reactive(), input$gene,
                                     style = input$gene_style %||% "raincloud",
                                     mode = mode) })
    })

    output$plot <- shiny::renderPlot({
      p <- cur_plot()
      if (inherits(p, "pheatmap")) { grid::grid.newpage(); grid::grid.draw(p$gtable) }
      else print(p)
    })

    output$qc_dl <- shiny::downloadHandler(
      filename = function() paste0("qc_", input$view %||% "plot", ".", input$qc_fmt),
      content  = function(file) {
        save_compare(cur_plot(), file, input$qc_fmt,
                     input$qc_w, input$qc_h, as.integer(input$qc_dpi))
      }
    )
  })
}
