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
                       " (genes × samples) and ",
                       shiny::strong("sample metadata"),
                       " (column 1 = sample ID matching counts colnames)."),
              shiny::p("Then run DESeq2 with your design variable."),
              shiny::p("Already have DE results from elsewhere? Skip to the third upload box."),
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
    bslib::nav_spacer(),
    bslib::nav_item(
      shiny::tags$a(
        href = "https://github.com/kmatmat/RNAflow",
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

  # DE analysis
  de_mod <- mod_de_server("de", data_mod)

  # Resolved DE results: prefer computed results, fall back to uploaded
  de_combined <- shiny::reactive({
    de_mod$de_results() %||% data_mod$de_results()
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
                "— using log2(counts+1) instead."),
          type = "warning", duration = 6
        )
        log2(counts + 1)
      }
    )
  })

  mod_volcano_server("volcano", de_combined)
  mod_heatmap_server("heatmap", de_combined, counts_norm, data_mod$metadata)
  mod_pca_server("pca", counts_norm, data_mod$metadata)
}
