#' Launch the RNAflow Shiny application
#'
#' Runs the full RNAflow Shiny app, which assembles all modules
#' (data, DE, volcano, heatmap, PCA, ...) into a single interface.
#'
#' @param port port to launch on (default: random free port)
#' @param launch_browser if TRUE, open in default browser
#' @return invisibly returns the Shiny app object
#' @export
#' @examples
#' \dontrun{
#'   RNAflow::run_app()
#' }
run_app <- function(port = NULL, launch_browser = TRUE) {
  shiny::shinyApp(
    ui = app_ui(),
    server = app_server,
    options = list(
      launch.browser = launch_browser,
      port = port
    )
  )
}

#' Alias for [run_app()]
#' @inheritParams run_app
#' @export
launch_app <- function(port = NULL, launch_browser = TRUE) {
  run_app(port = port, launch_browser = launch_browser)
}

#' RNAflow main UI
#'
#' @keywords internal
app_ui <- function() {
  bslib::page_navbar(
    title = shiny::tagList(
      shiny::span(class = "rnaflow-brand", "RNAflow"),
      shiny::tags$small(
        style = "font-weight:normal;color:#7F8C8D;margin-left:6px;",
        paste0("v", utils::packageVersion("RNAflow"))
      )
    ),
    theme = bslib::bs_theme(version = 5, bootswatch = "flatly",
                            primary = "#1D9E75"),
    header = shiny::tagList(
      shiny::tags$head(
        shiny::tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")
      )
    ),
    fillable = TRUE,
    bslib::nav_panel(
      title = "Data",
      shiny::fluidRow(
        shiny::column(4, mod_data_ui("data")),
        shiny::column(4, mod_de_ui("de")),
        shiny::column(
          4,
          bslib::card(
            bslib::card_header("Getting started"),
            bslib::card_body(
              shiny::p("Upload a ", shiny::strong("counts matrix"),
                       " (genes x samples) and ",
                       shiny::strong("sample metadata"),
                       " (column 1 = sample ID matching counts colnames)."),
              shiny::p("Then run DESeq2 with your design variable."),
              shiny::p("Already have DE results from elsewhere? Skip to the third upload box."),
              shiny::tags$hr(),
              shiny::p(shiny::tags$small(
                shiny::strong("Demo data"), " (in the package's ",
                shiny::tags$code("inst/extdata/"), " folder):",
                shiny::tags$ul(
                  shiny::tags$li(shiny::tags$code("demo_airway_counts.csv"),
                    " -- a real, published human dataset (airway, ",
                    "dexamethasone vs. control across 4 cell lines). Organism = ",
                    shiny::strong("Human"), "; design ",
                    shiny::tags$code("condition"), ", adjust for ",
                    shiny::tags$code("cell"), "."),
                  shiny::tags$li(shiny::tags$code("demo_multi_counts.csv"),
                    " -- a simulated mouse factorial set (6 groups) for ",
                    "multi-contrast / WGCNA testing. Organism = ",
                    shiny::strong("Mouse"), "; design ",
                    shiny::tags$code("group"), "."),
                  shiny::tags$li(shiny::tags$code("demo_counts.csv"),
                    " -- a minimal 2-group simulated set."))
              )),
              shiny::tags$hr(),
              shiny::p(shiny::tags$small(shiny::tags$em(
                "RNAflow validates everything on import. ",
                "If a file is malformed, you'll see exactly why."
              )))
            )
          )
        )
      )
    ),
    bslib::nav_panel("Volcano", mod_volcano_ui("volcano")),
    bslib::nav_panel("Heatmap", mod_heatmap_ui("heatmap")),
    bslib::nav_panel("PCA",     mod_pca_ui("pca")),
    bslib::nav_panel("QC",      mod_qc_ui("qc")),
    bslib::nav_panel("Compare", mod_compare_ui("compare")),
    bslib::nav_panel("Enrichment", mod_enrich_ui("enrich")),
    bslib::nav_panel("Network", mod_wgcna_ui("wgcna")),
    bslib::nav_panel("Activity", mod_activity_ui("activity")),
    bslib::nav_panel("AI", mod_ai_ui("ai")),
    bslib::nav_panel("Project", mod_project_ui("project")),
    bslib::nav_panel("Report", mod_report_ui("report")),
    bslib::nav_spacer(),
    bslib::nav_item(
      shiny::div(
        style = "display:flex;align-items:center;gap:6px;min-width:260px;",
        shiny::tags$small(style = "color:#7F8C8D;white-space:nowrap;", "Active:"),
        shiny::uiOutput("active_contrast_ui")
      )
    ),
    bslib::nav_item(
      shiny::tags$a(
        href = "https://github.com/KmBioChemo/RNAflow",
        target = "_blank",
        shiny::icon("github"), " GitHub"
      )
    )
  )
}

#' RNAflow main server
#'
#' @param input,output,session standard Shiny server args
#' @keywords internal
app_server <- function(input, output, session) {

  # Shared data layer
  data_mod <- mod_data_server("data")

  # Named store of DE contrasts, grown by each DESeq2 run
  contrasts_rv <- shiny::reactiveVal(list())

  # Enrichment / WGCNA settings captured for reproducible export
  settings_rv <- shiny::reactiveVal(list())

  # DE analysis (auto-adds each run to the store)
  de_mod <- mod_de_server("de", data_mod, contrasts_rv)

  # Mirror an uploaded pre-computed DE table into the store as well
  shiny::observeEvent(data_mod$de_results(), {
    res <- data_mod$de_results()
    if (is.null(res)) return()
    contrasts_rv(contrast_store_upsert(
      contrasts_rv(), "uploaded DE results", res,
      params = list(source = "upload")))
  })

  # Active-contrast selector shown in the navbar
  output$active_contrast_ui <- shiny::renderUI({
    store <- contrasts_rv()
    if (length(store) == 0) {
      return(shiny::tags$small(style = "color:#95A5A6;", "no contrast yet"))
    }
    shiny::selectInput("active_contrast", NULL,
                       choices = names(store),
                       selected = isolate_active(input$active_contrast, names(store)),
                       width = "210px")
  })

  # Resolved single DE table for the volcano / heatmap / PCA tabs:
  # the selected store entry, falling back to the latest computed/uploaded.
  de_combined <- shiny::reactive({
    store <- contrasts_rv()
    if (length(store) > 0) {
      sel <- isolate_active(input$active_contrast, names(store))
      return(store[[sel]]$results)
    }
    de_mod$de_results() %||% data_mod$de_results()
  })

  # Parameters of the active contrast (design_var / treated / reference), used
  # by the Heatmap / PCA "restrict to contrast" option.
  active_contrast_params <- shiny::reactive({
    store <- contrasts_rv()
    if (length(store) == 0) return(NULL)
    store[[isolate_active(input$active_contrast, names(store))]]$params
  })

  # Normalized counts for heatmap / PCA
  counts_norm <- shiny::reactive({
    counts <- data_mod$counts()
    if (is.null(counts)) return(NULL)
    tryCatch(
      normalize_counts(counts, data_mod$metadata(), method = "vst"),
      error = function(e) {
        shiny::showNotification(
          paste("Normalization failed:", conditionMessage(e),
                "-- using log2(counts+1) instead."),
          type = "warning", duration = 6
        )
        log2(counts + 1)
      }
    )
  })

  mod_volcano_server("volcano", de_combined)
  mod_heatmap_server("heatmap", de_combined, counts_norm, data_mod$metadata,
                     active_contrast_params)
  mod_pca_server("pca", counts_norm, data_mod$metadata, active_contrast_params)
  mod_qc_server("qc", de_combined, data_mod$counts, counts_norm, data_mod$metadata)
  mod_compare_server("compare", shiny::reactive(contrasts_rv()))
  enrich_result <- mod_enrich_server("enrich", de_combined, data_mod$organism,
                                     settings_rv)
  mod_wgcna_server("wgcna", counts_norm, data_mod$metadata, data_mod$organism,
                   settings_rv)
  mod_activity_server("activity", de_combined, data_mod$organism)
  mod_ai_server("ai", de_combined, enrich_result, data_mod$organism,
                active_contrast_params, settings_rv)
  mod_project_server("project", data_mod, contrasts_rv, settings_rv)
  mod_report_server("report", data_mod, contrasts_rv, settings_rv)
}

#' Resolve the active contrast label against the available choices
#'
#' Keeps the current selection if still valid, otherwise falls back to the
#' first available contrast. Avoids a transient NULL when the store changes.
#'
#' @param current the current `input$active_contrast` (may be NULL)
#' @param choices available contrast labels
#' @return a single valid label
#' @keywords internal
isolate_active <- function(current, choices) {
  if (length(choices) == 0) return(NULL)
  if (!is.null(current) && current %in% choices) return(current)
  choices[[1]]
}
