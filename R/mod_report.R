#' Reproducibility / report module
#'
#' Exports the session as a reproducible R script and a self-contained HTML
#' report, and shows the package versions used. Wraps [generate_r_script()]
#' and [build_report_html()].
#'
#' @param id namespace ID
#' @name mod_report
NULL

#' @rdname mod_report
#' @export
mod_report_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    bslib::layout_columns(
      col_widths = c(6, 6),
      bslib::card(
        bslib::card_header("Reproducible R script"),
        bslib::card_body(
          shiny::p(shiny::tags$small(
            "A runnable .R script that reproduces the whole pipeline ",
            "(DE, comparison, enrichment, network) with RNAflow's API -- ",
            "ready for a Methods section.")),
          shiny::downloadButton(ns("dl_script"), "Download .R script",
                                class = "btn btn-primary btn-sm")
        )
      ),
      bslib::card(
        bslib::card_header("HTML report"),
        bslib::card_body(
          shiny::p(shiny::tags$small(
            "A self-contained HTML report: parameters, DE summary, volcano ",
            "and cross-contrast figures, the reproducible script, and session ",
            "info. No external files -- open it anywhere.")),
          shiny::downloadButton(ns("dl_html"), "Download HTML report",
                                class = "btn btn-success btn-sm"),
          shiny::uiOutput(ns("report_note"))
        )
      )
    ),
    bslib::card(
      bslib::card_header("Script preview"),
      bslib::card_body(
        shiny::tags$pre(
          style = "max-height:340px;overflow:auto;font-size:12px;",
          shiny::textOutput(ns("script_preview"), container = shiny::tags$code)
        )
      )
    ),
    bslib::card(
      bslib::card_header("Session packages"),
      bslib::card_body(shiny::tableOutput(ns("manifest")))
    )
  )
}

#' @rdname mod_report
#' @param data_mod the value returned by [mod_data_server()]
#' @param contrast_store a `reactiveVal` holding the contrast store
#' @export
mod_report_server <- function(id, data_mod, contrast_store) {
  shiny::moduleServer(id, function(input, output, session) {

    project <- shiny::reactive({
      assemble_project(
        name      = "RNAflow analysis",
        organism  = data_mod$organism(),
        counts    = data_mod$counts(),
        metadata  = data_mod$metadata(),
        contrasts = contrast_store())
    })

    stamp <- function() format(Sys.time(), "%Y-%m-%d %H:%M")

    output$script_preview <- shiny::renderText({
      generate_r_script(project(), generated = stamp())
    })

    output$report_note <- shiny::renderUI({
      if (length(contrast_store()) == 0) {
        shiny::div(class = "demo-banner", style = "margin-top:8px;",
                   "Tip: run at least one contrast for a richer report.")
      }
    })

    output$manifest <- shiny::renderTable(session_manifest())

    output$dl_script <- shiny::downloadHandler(
      filename = function() "rnaflow_analysis.R",
      content  = function(file) {
        writeLines(generate_r_script(project(), generated = stamp()), file)
      }
    )

    output$dl_html <- shiny::downloadHandler(
      filename = function() "rnaflow_report.html",
      content  = function(file) {
        shiny::withProgress(message = "Building report...", value = 0.4, {
          build_report_html(project(), file, generated = stamp())
        })
      }
    )
  })
}
