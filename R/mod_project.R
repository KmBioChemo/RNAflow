#' Project manager module
#'
#' Save the full analysis session (counts, metadata, organism, and the whole
#' contrast store) to a `.rnaflow.rds` file, reload one, and re-open recent
#' projects. Restoration pushes state back into the data layer and the
#' contrast store so every downstream tab updates.
#'
#' @param id namespace ID
#' @name mod_project
NULL

#' @rdname mod_project
#' @export
mod_project_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_columns(
    col_widths = c(4, 4, 4),
    bslib::card(
      bslib::card_header("Save project"),
      bslib::card_body(
        shiny::textInput(ns("save_name"), "Project name", value = "my_project"),
        shiny::p(shiny::tags$small(shiny::tags$em(
          "Bundles counts, metadata, organism and all saved contrasts ",
          "into a single .rnaflow.rds file."))),
        shiny::downloadButton(ns("save_dl"), "Download project",
                              class = "btn btn-primary btn-sm"),
        shiny::uiOutput(ns("save_summary"))
      )
    ),
    bslib::card(
      bslib::card_header("Open project"),
      bslib::card_body(
        shiny::fileInput(ns("load_file"), "Load a .rnaflow.rds file",
                         accept = ".rds", buttonLabel = "Browse",
                         placeholder = "No file"),
        shiny::uiOutput(ns("load_status"))
      )
    ),
    bslib::card(
      bslib::card_header("Recent projects"),
      bslib::card_body(
        shiny::uiOutput(ns("recent_list"))
      )
    )
  )
}

#' @rdname mod_project
#' @param data_mod the value returned by [mod_data_server()] (must expose
#'   `set_state`)
#' @param contrast_store a `reactiveVal` holding the contrast store
#' @param settings_store optional `reactiveVal` with enrichment / WGCNA settings
#' @export
mod_project_server <- function(id, data_mod, contrast_store,
                               settings_store = NULL) {
  shiny::moduleServer(id, function(input, output, session) {

    recent_tick <- shiny::reactiveVal(0)  # bump to refresh the recent list

    gather_project <- function(name) {
      # Delegate to assemble_project() -- the single canonical constructor
      # (also used by the report module) so a saved project carries the FULL
      # session: enrichment, WGCNA, activity, signatures, AI interpretation and
      # the normalization method, not just enrichment + WGCNA.
      p <- assemble_project(
        name      = if (nzchar(name)) name else "untitled",
        organism  = data_mod$organism(),
        counts    = data_mod$counts(),
        metadata  = data_mod$metadata(),
        contrasts = contrast_store(),
        settings  = if (is.null(settings_store)) list() else settings_store())
      # Keep an uploaded DE table when no contrast was stored.
      if (is.null(p$de_results)) p$de_results <- data_mod$de_results()
      p
    }

    restore_project <- function(p) {
      # de_results kept NULL on purpose: the contrast store already holds
      # everything, and a non-NULL upload would spawn a duplicate entry.
      data_mod$set_state(counts = p$counts, metadata = p$metadata,
                         organism = p$organism, de_results = NULL)
      contrast_store(if (is.null(p$contrasts)) list() else p$contrasts)
      if (!is.null(settings_store)) {
        # Restore every settings slot so activity / signatures / AI / the
        # normalization method survive a save-reload, matching what was saved.
        settings_store(list(
          enrichment = p$enrichment %||% NULL,
          wgcna      = p$wgcna %||% NULL,
          activity   = p$activity %||% NULL,
          signatures = p$signatures %||% NULL,
          ai_interpretation = p$ai_interpretation %||% NULL,
          normalization_method = p$normalization_method %||% NULL))
      }
    }

    output$save_summary <- shiny::renderUI({
      store <- contrast_store()
      has_counts <- !is.null(data_mod$counts())
      shiny::div(class = "demo-banner", style = "margin-top:8px;",
                 sprintf("Counts: %s - %d contrast%s",
                         if (has_counts) "yes" else "no",
                         length(store), if (length(store) == 1) "" else "s"))
    })

    output$save_dl <- shiny::downloadHandler(
      filename = function() {
        nm <- gsub("[^A-Za-z0-9_-]+", "_", input$save_name %||% "project")
        paste0(if (nzchar(nm)) nm else "project", ".rnaflow.rds")
      },
      content = function(file) {
        p <- gather_project(input$save_name %||% "untitled")
        p$modified_at <- Sys.time()
        saveRDS(p, file)
        cache_recent_project(file, p$name)
        recent_tick(recent_tick() + 1)
      }
    )

    # Upload + restore
    shiny::observeEvent(input$load_file, {
      tryCatch({
        p <- load_project(input$load_file$datapath)
        restore_project(p)
        cache_recent_project(input$load_file$datapath, p$name)
        recent_tick(recent_tick() + 1)
        output$load_status <- shiny::renderUI(
          shiny::div(class = "demo-banner",
                     sprintf(" Loaded '%s' (%d contrasts)",
                             p$name, length(p$contrasts %||% list()))))
        shiny::showNotification(sprintf("Project '%s' loaded.", p$name),
                                type = "message", duration = 4)
      }, error = function(e) {
        output$load_status <- shiny::renderUI(
          shiny::div(class = "demo-banner",
                     style = "color:#C0392B;",
                     paste("Load failed:", conditionMessage(e))))
      })
    })

    # Recent projects list (rebuilt when recent_tick changes)
    output$recent_list <- shiny::renderUI({
      recent_tick()
      ns <- session$ns
      df <- list_recent_projects()
      if (nrow(df) == 0) {
        return(shiny::div(class = "demo-banner",
                          "No recent projects yet. Saved and opened projects ",
                          "appear here."))
      }
      items <- lapply(seq_len(nrow(df)), function(i) {
        token <- basename(df$file[i])
        onclick <- sprintf(
          "Shiny.setInputValue('%s', '%s', {priority:'event'});",
          ns("open_recent"), token)
        shiny::div(
          style = paste("display:flex;justify-content:space-between;",
                        "align-items:center;padding:5px 0;",
                        "border-bottom:1px solid #eee;"),
          shiny::div(
            shiny::strong(df$name[i]),
            shiny::tags$br(),
            shiny::tags$small(style = "color:#95A5A6;",
                              format(df$modified_at[i], "%Y-%m-%d %H:%M"))
          ),
          shiny::tags$button("Open", class = "btn btn-sm btn-outline-secondary",
                             onclick = onclick)
        )
      })
      shiny::div(items)
    })

    shiny::observeEvent(input$open_recent, {
      token <- basename(input$open_recent)  # guard against path traversal
      path <- file.path(rnaflow_recent_dir(), token)
      if (!file.exists(path)) {
        shiny::showNotification("That project file is no longer available.",
                                type = "warning")
        return()
      }
      tryCatch({
        p <- load_project(path)
        restore_project(p)
        shiny::showNotification(sprintf("Project '%s' opened.", p$name),
                                type = "message", duration = 4)
      }, error = function(e) {
        shiny::showNotification(paste("Open failed:", conditionMessage(e)),
                                type = "error", duration = 8)
      })
    })
  })
}
