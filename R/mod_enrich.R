#' Functional enrichment module
#'
#' Shiny module wrapping the [analysis_enrich] layer. Runs GSEA (against
#' MSigDB collections) or ORA (GO / KEGG / Reactome) on the active contrast,
#' and renders dotplot / bar / GSEA curve plus a results table.
#'
#' @param id namespace ID
#' @name mod_enrich
NULL

# Friendly collection name -> (collection, subcollection) for msigdbr.
ENRICH_COLLECTIONS <- list(
  "MSigDB Hallmark"             = list(collection = "H",  sub = NULL),
  "Reactome pathways (C2)"      = list(collection = "C2", sub = "CP:REACTOME"),
  "KEGG pathways (C2)"          = list(collection = "C2", sub = "CP:KEGG_LEGACY"),
  "GO Biological Process (C5)"  = list(collection = "C5", sub = "GO:BP"),
  "GO Molecular Function (C5)"  = list(collection = "C5", sub = "GO:MF"),
  "GO Cellular Component (C5)"  = list(collection = "C5", sub = "GO:CC")
)

#' @rdname mod_enrich
#' @export
mod_enrich_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 340,
      shiny::uiOutput(ns("organism_note")),
      shiny::radioButtons(ns("method"), "Method",
                          choices = c("GSEA" = "gsea", "ORA" = "ora"),
                          selected = "gsea", inline = TRUE),

      # ---- GSEA controls ----
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'gsea'", ns("method")),
        shiny::selectInput(ns("collection"), "Gene set collection",
                           choices = names(ENRICH_COLLECTIONS),
                           selected = "MSigDB Hallmark"),
        shiny::selectInput(ns("rank_by"), "Rank genes by",
                           choices = c("Wald statistic" = "stat",
                                       "signed -log10(p)" = "signed_p",
                                       "log2 fold change" = "log2fc"),
                           selected = "stat")
      ),

      # ---- ORA controls ----
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'ora'", ns("method")),
        shiny::selectInput(ns("db"), "Database",
                           choices = c("GO", "KEGG", "Reactome"),
                           selected = "GO"),
        shiny::conditionalPanel(
          sprintf("input['%s'] == 'GO'", ns("db")),
          shiny::selectInput(ns("ont"), "GO ontology",
                             choices = c("BP", "MF", "CC"), selected = "BP")
        ),
        shiny::radioButtons(ns("direction"), "Genes",
                            choices = c("Either" = "either", "Up" = "up",
                                        "Down" = "down"),
                            selected = "either", inline = TRUE),
        ui_slider_num(ns("lfc_sld"), ns("lfc_num"), "|log2FC| >", 0, 5, 1, 0.1),
        ui_slider_num(ns("padj_sld"), ns("padj_num"), "FDR (padj) <",
                      0.001, 0.2, 0.05, 0.001)
      ),

      shiny::tags$hr(style = "margin:8px 0;"),
      shiny::actionButton(ns("run"), "Run enrichment",
                          class = "btn btn-primary btn-sm",
                          style = "width:100%;"),
      shiny::tags$hr(style = "margin:8px 0;"),

      ui_section_title("Display"),
      shiny::selectInput(ns("view"), "Plot",
                         choices = c("Dotplot" = "dot", "Bar (-log10 FDR)" = "bar")),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'gsea' && input['%s'] == 'curve'",
                ns("method"), ns("view")),
        shiny::uiOutput(ns("pathway_pick"))
      ),
      ui_slider_num(ns("nterm_sld"), ns("nterm_num"), "Top terms", 5, 40, 20, 1),
      shiny::checkboxInput(ns("publication"), "Publication mode (export)", FALSE),
      ui_export_bar(ns("enr"), 7, 6)
    ),
    shiny::uiOutput(ns("notice")),
    shinycssloaders::withSpinner(
      shiny::plotOutput(ns("plot"), height = "520px"),
      type = 6, color = "#1D9E75"
    ),
    shiny::tags$hr(),
    DT::DTOutput(ns("table"))
  )
}

#' @rdname mod_enrich
#' @param de_reactive a reactive returning the active contrast DE data.frame
#' @param organism_reactive a reactive returning the organism keyword
#' @export
mod_enrich_server <- function(id, de_reactive, organism_reactive) {
  shiny::moduleServer(id, function(input, output, session) {

    mirror <- function(s, n) {
      shiny::observeEvent(input[[s]],
        shiny::updateNumericInput(session, n, value = input[[s]]), ignoreInit = TRUE)
      shiny::observeEvent(input[[n]],
        shiny::updateSliderInput(session, s, value = input[[n]]), ignoreInit = TRUE)
    }
    mirror("lfc_sld", "lfc_num"); mirror("padj_sld", "padj_num")
    mirror("nterm_sld", "nterm_num")

    result <- shiny::reactiveVal(NULL)   # list(method, table, gene_sets, de)

    output$organism_note <- shiny::renderUI({
      org <- organism_reactive() %||% "human"
      shiny::div(class = "demo-banner", style = "margin-bottom:6px;",
                 sprintf("Organism: %s (set on the Data tab)", org))
    })

    # Offer the GSEA curve view only after a GSEA run
    shiny::observeEvent(result(), {
      r <- result()
      choices <- c("Dotplot" = "dot", "Bar (-log10 FDR)" = "bar")
      if (!is.null(r) && r$method == "gsea") {
        choices <- c(choices, "GSEA curve" = "curve")
      }
      shiny::updateSelectInput(session, "view", choices = choices,
                               selected = shiny::isolate(input$view))
    })

    shiny::observeEvent(input$run, {
      de  <- de_reactive()
      org <- organism_reactive() %||% "human"
      if (is.null(de)) {
        shiny::showNotification("No active contrast. Run DESeq2 first.",
                                type = "warning"); return()
      }
      shiny::withProgress(message = "Running enrichment...", value = 0.3, {
        tryCatch({
          if (input$method == "gsea") {
            cc <- ENRICH_COLLECTIONS[[input$collection]]
            gs <- get_gene_sets(org, collection = cc$collection,
                                subcollection = cc$sub, id_type = "symbol")
            shiny::incProgress(0.4)
            tab <- run_gsea(de, gs, rank_by = input$rank_by, min_size = 15)
            result(list(method = "gsea", table = tab, gene_sets = gs, de = de))
          } else {
            sig <- contrast_sig_genes(de, input$padj_num %||% 0.05,
                                      input$lfc_num %||% 1, input$direction)
            if (length(sig) < 5) stop("Only ", length(sig),
                                      " significant genes — loosen thresholds.")
            shiny::incProgress(0.4)
            tab <- run_ora(sig, org, db = input$db,
                           ont = if (input$db == "GO") input$ont else "BP",
                           universe = de$gene)
            result(list(method = "ora", table = tab, gene_sets = NULL, de = de))
          }
          shiny::showNotification(
            sprintf("Enrichment done: %d term%s.", nrow(result()$table),
                    if (nrow(result()$table) == 1) "" else "s"),
            type = "message", duration = 4)
        }, error = function(e) {
          result(NULL)
          shiny::showNotification(paste("Enrichment failed:", conditionMessage(e)),
                                  type = "error", duration = 10)
        })
      })
    })

    output$pathway_pick <- shiny::renderUI({
      r <- result()
      shiny::req(r, r$method == "gsea", nrow(r$table) > 0)
      shiny::selectInput(session$ns("pathway"), "Pathway (curve)",
                         choices = r$table$pathway,
                         selected = r$table$pathway[1])
    })

    output$notice <- shiny::renderUI({
      r <- result()
      if (is.null(r)) {
        return(shiny::div(class = "demo-banner", style = "margin-bottom:8px;",
                          "Pick a method and click Run enrichment."))
      }
      if (nrow(r$table) == 0) {
        return(shiny::div(class = "demo-banner", style = "margin-bottom:8px;",
                          "No significant terms at the current settings."))
      }
      NULL
    })

    cur_plot <- shiny::reactive({
      r <- result()
      shiny::req(r, nrow(r$table) > 0)
      mode <- if (isTRUE(input$publication)) "publication" else "exploration"
      n <- max(5L, as.integer(input$nterm_num %||% 20))
      if (r$method == "gsea" && (input$view %||% "dot") == "curve") {
        shiny::req(input$pathway)
        return(fig_gsea_curve(r$de, r$gene_sets[[input$pathway]],
                              rank_by = input$rank_by, title = input$pathway,
                              mode = mode))
      }
      if ((input$view %||% "dot") == "bar") fig_enrich_bar(r$table, n = n, mode = mode)
      else fig_enrich_dot(r$table, n = n, mode = mode)
    })

    output$plot <- shiny::renderPlot({ print(cur_plot()) })

    output$table <- DT::renderDT({
      r <- result()
      shiny::req(r, nrow(r$table) > 0)
      df <- r$table
      if (r$method == "gsea") {
        df <- df[, intersect(c("pathway", "NES", "size", "pval", "padj",
                               "leading_edge"), colnames(df))]
        df$NES <- round(df$NES, 3)
      } else {
        df <- df[, intersect(c("Description", "Count", "GeneRatio",
                               "pvalue", "padj"), colnames(df))]
      }
      for (col in c("pval", "pvalue", "padj")) {
        if (col %in% colnames(df)) df[[col]] <- formatC(df[[col]], format = "e", digits = 2)
      }
      DT::datatable(df, rownames = FALSE, class = "compact stripe hover",
                    options = list(pageLength = 10, scrollX = TRUE))
    })

    output$enr_dl <- shiny::downloadHandler(
      filename = function() paste0("enrichment_", result()$method %||% "plot",
                                   ".", input$enr_fmt),
      content = function(file) {
        save_ggplot(cur_plot(), file, input$enr_fmt,
                    input$enr_w, input$enr_h, as.integer(input$enr_dpi))
      }
    )
  })
}
