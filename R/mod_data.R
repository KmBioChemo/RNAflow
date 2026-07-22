#' Data input module
#'
#' Handles file upload (counts, metadata, pre-computed DE results) and
#' exposes validated reactive objects to other modules.
#'
#' Returns a list of reactives: `counts()`, `metadata()`, `de_results()`,
#' `organism()`.
#'
#' @param id namespace ID
#' @name mod_data
NULL

#' @rdname mod_data
#' @export
mod_data_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::card(
    bslib::card_header("Data input"),
    bslib::card_body(
      ui_section_title("Try a demo dataset"),
      shiny::div(
        style = "display:flex; flex-direction:column; gap:6px;",
        shiny::actionButton(
          ns("demo_airway"), "Airway (8 samples)",
          icon = shiny::icon("flask"),
          class = "btn btn-outline-primary btn-sm",
          width = "100%"),
        shiny::actionButton(
          ns("demo_tcga"), "TCGA pan-cancer (120 samples)",
          icon = shiny::icon("dna"),
          class = "btn btn-outline-primary btn-sm",
          width = "100%")
      ),
      shiny::p(class = "rf-microcopy",
               "Real published data bundled with the package — loads counts ",
               "and metadata in one click, no download needed."),
      shiny::tags$hr(style = "margin:8px 0;"),
      shiny::selectInput(ns("organism"), "Organism",
                         choices = c("Human" = "human",
                                     "Mouse" = "mouse",
                                     "Rat"   = "rat"),
                         selected = "human"),
      shiny::tags$hr(style = "margin:8px 0;"),
      ui_section_title("Counts matrix (genes x samples)"),
      shiny::fileInput(ns("counts_file"), NULL,
                       accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls"),
                       buttonLabel = "Browse", placeholder = "No file"),
      shiny::uiOutput(ns("counts_status")),
      shiny::tags$hr(style = "margin:8px 0;"),
      ui_section_title("Sample metadata"),
      shiny::fileInput(ns("meta_file"), NULL,
                       accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls"),
                       buttonLabel = "Browse", placeholder = "No file"),
      shiny::uiOutput(ns("meta_status")),
      shiny::tags$hr(style = "margin:8px 0;"),
      ui_section_title("OR upload pre-computed DE results"),
      shiny::fileInput(ns("de_file"), NULL,
                       accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls"),
                       buttonLabel = "Browse", placeholder = "No file"),
      shiny::uiOutput(ns("de_status"))
    )
  )
}

#' @rdname mod_data
#' @export
mod_data_server <- function(id) {
  shiny::moduleServer(id, function(input, output, session) {

    counts_r <- shiny::reactiveVal(NULL)
    meta_r   <- shiny::reactiveVal(NULL)
    de_r     <- shiny::reactiveVal(NULL)

    # Bundled demo datasets: load the packaged counts + metadata CSVs through
    # the same validated readers as an upload. Files live in inst/extdata/ and
    # are resolved with system.file() so they work from an installed package
    # (there is no browsable inst/extdata/ folder after install).
    load_demo <- function(base, organism) {
      cf <- system.file("extdata", paste0(base, "_counts.csv"),
                        package = "RNAflow")
      mf <- system.file("extdata", paste0(base, "_metadata.csv"),
                        package = "RNAflow")
      if (!nzchar(cf) || !nzchar(mf)) {
        shiny::showNotification(
          "Demo files were not found in the installed package.",
          type = "error", duration = 8)
        return(invisible(NULL))
      }
      tryCatch({
        m <- read_counts(cf, ext = "csv",
                         validate = TRUE, strict_integer = TRUE)
        d <- read_metadata(mf, ext = "csv", validate = TRUE,
                           counts_samples = colnames(m))
        counts_r(m)
        meta_r(d)
        de_r(NULL)
        shiny::updateSelectInput(session, "organism", selected = organism)
        shiny::showNotification(
          sprintf("Demo loaded: %d genes × %d samples.",
                  nrow(m), ncol(m)),
          type = "message", duration = 5)
      }, error = function(e) {
        shiny::showNotification(
          paste("Demo load failed:", conditionMessage(e)),
          type = "error", duration = 10)
      })
    }

    shiny::observeEvent(input$demo_airway, load_demo("demo_airway", "human"))
    shiny::observeEvent(input$demo_tcga,   load_demo("demo_tcga", "human"))

    # Counts upload
    shiny::observeEvent(input$counts_file, {
      tryCatch({
        m <- read_counts(input$counts_file$datapath,
                         ext = tolower(tools::file_ext(input$counts_file$name)),
                         validate = TRUE, strict_integer = TRUE)
        counts_r(m)
        shiny::showNotification(
          sprintf("Counts loaded: %d genes x %d samples",
                  nrow(m), ncol(m)),
          type = "message", duration = 4
        )
      }, error = function(e) {
        counts_r(NULL)
        shiny::showNotification(paste("Counts upload failed:", conditionMessage(e)),
                                type = "error", duration = 10)
      })
    })

    # Metadata upload
    shiny::observeEvent(input$meta_file, {
      tryCatch({
        df <- read_metadata(
          input$meta_file$datapath,
          ext = tolower(tools::file_ext(input$meta_file$name)),
          validate = TRUE,
          counts_samples = if (!is.null(counts_r())) colnames(counts_r()) else NULL
        )
        meta_r(df)
        shiny::showNotification(
          sprintf("Metadata loaded: %d samples x %d annotations",
                  nrow(df), ncol(df) - 1),
          type = "message", duration = 4
        )
      }, error = function(e) {
        meta_r(NULL)
        shiny::showNotification(paste("Metadata upload failed:", conditionMessage(e)),
                                type = "error", duration = 10)
      })
    })

    # Pre-computed DE results
    shiny::observeEvent(input$de_file, {
      tryCatch({
        df <- read_de_results(
          input$de_file$datapath,
          ext = tolower(tools::file_ext(input$de_file$name))
        )
        de_r(df)
        shiny::showNotification(
          sprintf("DE results loaded: %d genes", nrow(df)),
          type = "message", duration = 4
        )
      }, error = function(e) {
        de_r(NULL)
        shiny::showNotification(paste("DE results upload failed:", conditionMessage(e)),
                                type = "error", duration = 10)
      })
    })

    # Status messages
    output$counts_status <- shiny::renderUI({
      m <- counts_r()
      if (is.null(m)) return(shiny::div(class = "demo-banner", "No counts loaded yet."))
      shiny::div(class = "demo-banner",
                 sprintf("\u2713 %d genes \u00D7 %d samples", nrow(m), ncol(m)))
    })
    output$meta_status <- shiny::renderUI({
      d <- meta_r()
      if (is.null(d)) return(shiny::div(class = "demo-banner", "No metadata loaded yet."))
      shiny::div(class = "demo-banner",
                 sprintf("\u2713 %d samples, columns: %s",
                         nrow(d), paste(colnames(d)[-1], collapse = ", ")))
    })
    output$de_status <- shiny::renderUI({
      d <- de_r()
      if (is.null(d)) return(NULL)
      shiny::div(class = "demo-banner",
                 sprintf("\u2713 DE results: %d genes", nrow(d)))
    })

    # Restore data layer from a loaded project (used by the project manager)
    set_state <- function(counts = NULL, metadata = NULL, organism = NULL,
                          de_results = NULL) {
      counts_r(counts)
      meta_r(metadata)
      de_r(de_results)
      if (!is.null(organism)) {
        shiny::updateSelectInput(session, "organism", selected = organism)
      }
    }

    list(
      counts     = shiny::reactive(counts_r()),
      metadata   = shiny::reactive(meta_r()),
      de_results = shiny::reactive(de_r()),
      organism   = shiny::reactive(input$organism),
      set_state  = set_state
    )
  })
}
