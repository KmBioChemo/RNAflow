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
      shiny::selectInput(ns("organism"), "Organism",
                         choices = c("Human" = "human",
                                     "Mouse" = "mouse",
                                     "Rat"   = "rat"),
                         selected = "human"),
      shiny::tags$hr(style = "margin:8px 0;"),
      ui_section_title("Counts matrix (genes × samples)"),
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

    # Counts upload
    shiny::observeEvent(input$counts_file, {
      tryCatch({
        m <- read_counts(input$counts_file$datapath,
                         ext = tolower(tools::file_ext(input$counts_file$name)),
                         validate = TRUE, strict_integer = TRUE)
        counts_r(m)
        shiny::showNotification(
          sprintf("Counts loaded: %d genes × %d samples",
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
          sprintf("Metadata loaded: %d samples × %d annotations",
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
