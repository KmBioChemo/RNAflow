#' Signatures module (per-sample GSVA / ssGSEA scores)
#'
#' A per-sample gene-set scoring tab: pick an MSigDB collection and a method
#' (GSVA or ssGSEA), score every sample, and view the sets x samples signature
#' matrix as an annotated heatmap. Thin Shiny wrapper over [run_gsva()] /
#' [fig_gsva_heatmap()] (both guarded on the optional \pkg{GSVA} dependency).
#'
#' @name mod_signatures
#' @keywords internal
NULL

#' @rdname mod_signatures
#' @param id module id
#' @export
mod_signatures_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    ui_page_header("Signatures (GSVA / ssGSEA)",
                   "Per-sample gene-set scores from normalized counts."),
    sidebar = bslib::sidebar(
      width = 300,
      shiny::selectInput(ns("collection"), "Gene set collection",
                         choices = names(ENRICH_COLLECTIONS),
                         selected = names(ENRICH_COLLECTIONS)[1]),
      shiny::radioButtons(ns("method"), "Method",
                          choices = c("GSVA" = "gsva", "ssGSEA" = "ssgsea"),
                          selected = "gsva", inline = TRUE),
      ui_slider_num(ns("nset_sld"), ns("nset_num"), "Top sets shown",
                    5, 60, 30, 5),
      shiny::uiOutput(ns("group_ui")),
      shiny::actionButton(ns("run"), "Compute scores",
                          class = "btn btn-primary",
                          icon = shiny::icon("play")),
      ui_export_bar(ns("gsva"), 8, 7)
    ),
    shiny::uiOutput(ns("notice")),
    shinycssloaders::withSpinner(
      shiny::plotOutput(ns("plot"), height = "600px"),
      type = 6, color = "#1D9E75")
  )
}

#' @rdname mod_signatures
#' @param counts_norm_reactive reactive: normalized counts matrix
#' @param metadata_reactive reactive: sample metadata
#' @param organism_reactive reactive: organism keyword
#' @param settings_store optional `reactiveVal`; the last run is recorded under
#'   `$gsva` for reproducibility
#' @export
mod_signatures_server <- function(id, counts_norm_reactive, metadata_reactive,
                                  organism_reactive, settings_store = NULL) {
  shiny::moduleServer(id, function(input, output, session) {

    shiny::observeEvent(input$nset_sld,
      shiny::updateNumericInput(session, "nset_num", value = input$nset_sld),
      ignoreInit = TRUE)
    shiny::observeEvent(input$nset_num,
      shiny::updateSliderInput(session, "nset_sld", value = input$nset_num),
      ignoreInit = TRUE)

    scores <- shiny::reactiveVal(NULL)

    output$group_ui <- shiny::renderUI({
      md <- metadata_reactive()
      ch <- if (!is.null(md) && ncol(md) >= 2) c("(none)", colnames(md)[-1]) else "(none)"
      shiny::selectInput(session$ns("group_by"), "Annotate samples by",
                         choices = ch,
                         selected = if (length(ch) > 1) ch[2] else "(none)")
    })

    have_deps <- function() {
      requireNamespace("GSVA", quietly = TRUE) &&
        requireNamespace("msigdbr", quietly = TRUE)
    }

    output$notice <- shiny::renderUI({
      if (!have_deps()) {
        return(ui_banner(
          "Per-sample scoring needs the 'GSVA' and 'msigdbr' packages. ",
          "Install with BiocManager::install(c('GSVA','msigdbr')).",
          type = "warning"))
      }
      if (is.null(counts_norm_reactive())) {
        return(ui_banner("Load a counts matrix (Data tab) to compute ",
                         "per-sample signatures.", type = "warning"))
      }
      if (is.null(scores())) {
        return(shiny::div(class = "demo-banner",
                          "Pick a collection and click ",
                          shiny::strong("Compute scores"), "."))
      }
      NULL
    })

    shiny::observeEvent(input$run, {
      cn <- counts_norm_reactive()
      if (is.null(cn)) {
        shiny::showNotification("No counts matrix loaded.", type = "warning")
        return()
      }
      if (!have_deps()) {
        shiny::showNotification("Install GSVA + msigdbr first.", type = "error")
        return()
      }
      org <- organism_reactive() %||% "human"
      cc  <- ENRICH_COLLECTIONS[[input$collection]]
      shiny::withProgress(message = "Scoring samples...", value = 0.3, {
        tryCatch({
          sets <- get_gene_sets(org, cc$collection, cc$sub)
          cm   <- gsva_symbol_counts(cn, org)
          es   <- run_gsva(cm, sets, method = input$method %||% "gsva")
          scores(es)
          if (!is.null(settings_store)) {
            s <- settings_store()
            s$gsva <- list(collection = input$collection,
                           method = input$method %||% "gsva",
                           n_sets = nrow(es), n_samples = ncol(es),
                           generated = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
            settings_store(s)
          }
        }, error = function(e) {
          scores(NULL)
          shiny::showNotification(paste("Scoring failed:", conditionMessage(e)),
                                  type = "error", duration = 10)
        })
      })
    })

    cur_hm <- shiny::reactive({
      es <- scores(); shiny::req(es)
      grp <- if (!is.null(input$group_by) && input$group_by != "(none)")
        input$group_by else NULL
      fig_gsva_heatmap(es, metadata_reactive(), group_by = grp,
                       n_top = as.integer(input$nset_num %||% 30))
    })

    output$plot <- shiny::renderPlot({
      grid::grid.newpage(); grid::grid.draw(cur_hm()$gtable)
    })

    output$gsva_dl <- shiny::downloadHandler(
      filename = function() paste0("gsva_scores.", input$gsva_fmt),
      content = function(file) {
        save_pheatmap(cur_hm(), file, input$gsva_fmt,
                      input$gsva_w, input$gsva_h, as.integer(input$gsva_dpi))
      })
  })
}
