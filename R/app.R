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
                            primary = "#1D9E75", success = "#1D9E75"),
    header = shiny::tagList(
      shiny::tags$head(
        shiny::tags$link(rel = "stylesheet", type = "text/css",
                         href = "rnaflow/rnaflow.css")
      )
    ),
    fillable = TRUE,
    bslib::nav_panel(
      title = "Data",
      icon = shiny::icon("table"),
      shiny::fluidRow(
        shiny::column(4, mod_data_ui("data")),
        shiny::column(4, mod_de_ui("de")),
        shiny::column(
          4,
          bslib::card(
            bslib::card_header("Getting started"),
            bslib::card_body(
              ui_section_title("Workflow"),
              shiny::tags$ol(
                class = "rf-steps",
                shiny::tags$li("Upload a ", shiny::strong("counts matrix"),
                               " (genes x samples) and ",
                               shiny::strong("sample metadata"),
                               " (column 1 = sample ID, matching the counts columns)."),
                shiny::tags$li("Run ", shiny::strong("DESeq2"),
                               " with your design variable."),
                shiny::tags$li("Explore, enrich, and export -- or load a saved project.")
              ),
              shiny::p(class = "rf-microcopy",
                       "Already have DE results from elsewhere? Skip straight to the ",
                       "third upload box."),
              ui_section_title("Demo datasets"),
              shiny::p(class = "rf-microcopy",
                       "Bundled in the package's ",
                       shiny::tags$code("inst/extdata/"), " folder:"),
              shiny::tags$ul(
                class = "rf-demo-list",
                shiny::tags$li(shiny::tags$code("demo_airway_counts.csv"),
                  " -- real published human data (Himes et al. 2014, airway ",
                  "smooth muscle; dexamethasone vs. control across 4 cell ",
                  "lines). Organism ", shiny::strong("Human"), "; design ",
                  shiny::tags$code("condition"), ", adjust for ",
                  shiny::tags$code("cell"), ". Gene symbols."),
                shiny::tags$li(shiny::tags$code("demo_pickrell_counts.csv"),
                  " -- real published human data (Pickrell et al. 2010, ",
                  "lymphoblastoid cell lines; a balanced female-vs-male subset). ",
                  "Organism ", shiny::strong("Human"), "; design ",
                  shiny::tags$code("sex"), ". Ensembl IDs.")
              ),
              shiny::div(
                class = "rnaflow-banner rf-success",
                shiny::icon("circle-check", class = "rf-ic"),
                shiny::span("RNAflow validates every file on import -- if ",
                            "something is malformed, you'll see exactly why.")
              )
            )
          )
        )
      )
    ),
    bslib::nav_panel("Volcano", icon = shiny::icon("chart-simple"),
                     mod_volcano_ui("volcano")),
    bslib::nav_panel("Explore", icon = shiny::icon("table-cells"),
                     mod_linked_ui("linked")),
    bslib::nav_panel("Heatmap", icon = shiny::icon("table-cells-large"),
                     mod_heatmap_ui("heatmap")),
    bslib::nav_panel("PCA", icon = shiny::icon("braille"),
                     mod_pca_ui("pca")),
    bslib::nav_panel("QC", icon = shiny::icon("magnifying-glass-chart"),
                     mod_qc_ui("qc")),
    bslib::nav_panel("Compare", icon = shiny::icon("code-compare"),
                     mod_compare_ui("compare")),
    bslib::nav_panel("Enrichment", icon = shiny::icon("diagram-project"),
                     mod_enrich_ui("enrich")),
    bslib::nav_panel("Network", icon = shiny::icon("share-nodes"),
                     mod_wgcna_ui("wgcna")),
    bslib::nav_panel("Activity", icon = shiny::icon("wave-square"),
                     mod_activity_ui("activity")),
    bslib::nav_panel("Signatures", icon = shiny::icon("fingerprint"),
                     mod_signatures_ui("signatures")),
    bslib::nav_panel("AI", icon = shiny::icon("robot"),
                     mod_ai_ui("ai")),
    bslib::nav_panel("Project", icon = shiny::icon("folder-open"),
                     mod_project_ui("project")),
    bslib::nav_panel("Report", icon = shiny::icon("file-lines"),
                     mod_report_ui("report")),
    bslib::nav_spacer(),
    bslib::nav_item(
      shiny::div(
        class = "rnaflow-active",
        shiny::span(class = "rf-active-label", "Active contrast"),
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
  mod_linked_server("linked", de_combined)
  mod_heatmap_server("heatmap", de_combined, counts_norm, data_mod$metadata,
                     active_contrast_params)
  mod_pca_server("pca", counts_norm, data_mod$metadata, active_contrast_params)
  mod_qc_server("qc", de_combined, data_mod$counts, counts_norm, data_mod$metadata)
  mod_compare_server("compare", shiny::reactive(contrasts_rv()))
  enrich_result <- mod_enrich_server("enrich", de_combined, data_mod$organism,
                                     settings_rv)
  mod_wgcna_server("wgcna", counts_norm, data_mod$metadata, data_mod$organism,
                   settings_rv)
  mod_activity_server("activity", de_combined, data_mod$organism, settings_rv)
  mod_signatures_server("signatures", counts_norm, data_mod$metadata,
                        data_mod$organism, settings_rv)
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
