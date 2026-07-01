#' WGCNA co-expression network module
#'
#' Shiny module wrapping the [analysis_wgcna] layer: soft-threshold picking,
#' module detection, module-trait correlation, hub genes, eigengene profiles,
#' and per-module pathway enrichment (reusing [run_ora()] from phase 3).
#'
#' @param id namespace ID
#' @name mod_wgcna
NULL

#' @rdname mod_wgcna
#' @export
mod_wgcna_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 340,
      ui_section_title("1. Network construction"),
      ui_slider_num(ns("ngenes_sld"), ns("ngenes_num"),
                    "Top variable genes", 500, 8000, 3000, 250),
      shiny::selectInput(ns("net_type"), "Network type",
                         choices = c("signed", "unsigned", "signed hybrid"),
                         selected = "signed"),
      shiny::actionButton(ns("pick"), "Pick soft threshold",
                          class = "btn btn-outline-primary btn-sm",
                          style = "width:100%;"),
      shiny::tags$hr(style = "margin:8px 0;"),
      ui_section_title("2. Module detection"),
      shiny::numericInput(ns("power"), "Soft-thresholding power",
                          value = 6, min = 1, max = 30, step = 1),
      ui_advanced_panel(
        shiny::numericInput(ns("min_size"), "Min module size", 30, 5, 200, 5),
        shiny::numericInput(ns("merge_cut"), "Merge cut height", 0.25, 0, 1, 0.05),
        shiny::numericInput(ns("deep_split"), "Deep split (0-4)", 2, 0, 4, 1)
      ),
      shiny::actionButton(ns("detect"), "Detect modules",
                          class = "btn btn-primary btn-sm", style = "width:100%;"),
      shiny::p(shiny::tags$small(
        style = "color:#7F8C8D;display:block;margin-top:6px;",
        "WGCNA is exploratory and most reliable with larger sample sizes ",
        "(roughly 15 or more). Treat modules from small datasets as ",
        "hypotheses.")),
      shiny::tags$hr(style = "margin:8px 0;"),
      ui_section_title("Display"),
      shiny::radioButtons(
        ns("view"), NULL,
        choices = c("Soft threshold" = "sft", "Module-trait heatmap" = "mt",
                    "Module sizes" = "sizes", "Module eigengene" = "eigen",
                    "Module network" = "network"),
        selected = "sft"),
      shiny::conditionalPanel(
        sprintf("input['%s'] == 'eigen' || input['%s'] == 'network'",
                ns("view"), ns("view")),
        shiny::uiOutput(ns("module_pick"))
      ),
      ui_export_bar(ns("net"), 8, 5)
    ),
    shiny::uiOutput(ns("notice")),
    shinycssloaders::withSpinner(
      shiny::plotOutput(ns("plot"), height = "480px"),
      type = 6, color = "#1D9E75"),
    shiny::conditionalPanel(
      sprintf("output['%s']", ns("has_modules")),
      shiny::tags$hr(),
      bslib::layout_columns(
        col_widths = c(6, 6),
        bslib::card(
          bslib::card_header("Hub genes (top by module membership)"),
          DT::DTOutput(ns("hub_table"))
        ),
        bslib::card(
          bslib::card_header(
            shiny::div(style = "display:flex;justify-content:space-between;align-items:center;gap:6px;",
                       shiny::span("Module enrichment (GO BP)"),
                       shiny::div(
                         shiny::actionButton(ns("enrich"), "This module",
                                             class = "btn btn-success btn-sm"),
                         shiny::actionButton(ns("enrich_all"), "All modules",
                                             class = "btn btn-outline-success btn-sm")))),
          shinycssloaders::withSpinner(
            shiny::plotOutput(ns("enrich_plot"), height = "340px"),
            type = 6, color = "#1D9E75"),
          DT::DTOutput(ns("enrich_table"))
        )
      )
    )
  )
}

#' @rdname mod_wgcna
#' @param counts_norm_reactive reactive returning the normalized counts matrix
#'   (genes x samples)
#' @param metadata_reactive reactive returning the sample metadata
#' @param organism_reactive reactive returning the organism keyword
#' @export
mod_wgcna_server <- function(id, counts_norm_reactive, metadata_reactive,
                             organism_reactive) {
  shiny::moduleServer(id, function(input, output, session) {

    shiny::observeEvent(input$ngenes_sld,
      shiny::updateNumericInput(session, "ngenes_num", value = input$ngenes_sld),
      ignoreInit = TRUE)
    shiny::observeEvent(input$ngenes_num,
      shiny::updateSliderInput(session, "ngenes_sld", value = input$ngenes_num),
      ignoreInit = TRUE)

    sft_rv     <- shiny::reactiveVal(NULL)
    wg_rv      <- shiny::reactiveVal(NULL)
    enrich_rv  <- shiny::reactiveVal(NULL)   # single-module ORA
    modules_rv <- shiny::reactiveVal(NULL)   # all-module comparison
    enrich_view <- shiny::reactiveVal("single")

    datexpr <- function() {
      cn <- counts_norm_reactive()
      shiny::validate(shiny::need(!is.null(cn),
        "Load counts and run DESeq2 first so a normalized matrix is available."))
      wgcna_datexpr(cn, n_genes = as.integer(input$ngenes_num %||% 3000))
    }

    # Step 1: soft threshold
    shiny::observeEvent(input$pick, {
      shiny::withProgress(message = "Picking soft threshold...", value = 0.4, {
        tryCatch({
          sft <- wgcna_pick_power(datexpr(), network_type = input$net_type)
          sft_rv(sft)
          shiny::updateNumericInput(session, "power", value = sft$suggested)
          shiny::updateRadioButtons(session, "view", selected = "sft")
          shiny::showNotification(
            sprintf("Suggested power: %s", sft$suggested),
            type = "message", duration = 4)
        }, error = function(e)
          shiny::showNotification(paste("Soft threshold failed:",
                                        conditionMessage(e)),
                                  type = "error", duration = 10))
      })
    })

    # Step 2: module detection
    shiny::observeEvent(input$detect, {
      shiny::withProgress(message = "Detecting modules...", value = 0.4, {
        tryCatch({
          wg <- run_wgcna(datexpr(), power = as.integer(input$power %||% 6),
                          network_type = input$net_type,
                          min_module_size = as.integer(input$min_size %||% 30),
                          merge_cut_height = as.numeric(input$merge_cut %||% 0.25),
                          deep_split = as.integer(input$deep_split %||% 2))
          wg_rv(wg); enrich_rv(NULL)
          shiny::updateRadioButtons(session, "view", selected = "mt")
          shiny::showNotification(
            sprintf("Found %d modules (+ grey).",
                    length(setdiff(unique(wg$modules), "grey"))),
            type = "message", duration = 4)
        }, error = function(e)
          shiny::showNotification(paste("Module detection failed:",
                                        conditionMessage(e)),
                                  type = "error", duration = 12))
      })
    })

    output$has_modules <- shiny::reactive(!is.null(wg_rv()))
    shiny::outputOptions(output, "has_modules", suspendWhenHidden = FALSE)

    output$module_pick <- shiny::renderUI({
      wg <- wg_rv(); shiny::req(wg)
      mods <- setdiff(unique(wg$modules), "grey")
      shiny::selectInput(session$ns("module"), "Module",
                         choices = mods, selected = mods[1])
    })

    traits_cor <- shiny::reactive({
      wg <- wg_rv(); md <- metadata_reactive()
      shiny::req(wg, md)
      tr <- build_traits(md, rownames(wg$MEs))
      module_trait_cor(wg$MEs, tr)
    })

    groups_vec <- shiny::reactive({
      md <- metadata_reactive(); wg <- wg_rv()
      shiny::req(md, wg)
      stats::setNames(as.character(md[[2]])[match(rownames(wg$MEs), md[[1]])],
                      rownames(wg$MEs))
    })

    output$notice <- shiny::renderUI({
      v <- input$view %||% "sft"
      if (v == "sft" && is.null(sft_rv())) {
        return(shiny::div(class = "demo-banner", style = "margin-bottom:8px;",
          "Step 1: choose the number of genes, then 'Pick soft threshold'."))
      }
      if (v != "sft" && is.null(wg_rv())) {
        return(shiny::div(class = "demo-banner", style = "margin-bottom:8px;",
          "Step 2: set the power and click 'Detect modules'."))
      }
      NULL
    })

    cur_plot <- shiny::reactive({
      v <- input$view %||% "sft"
      mode <- "exploration"
      if (v == "sft") { shiny::req(sft_rv()); return(fig_soft_threshold(sft_rv(), mode)) }
      wg <- wg_rv(); shiny::req(wg)
      switch(v,
        mt      = fig_module_trait(traits_cor(), mode),
        sizes   = fig_module_sizes(wg, mode = mode),
        eigen   = { shiny::req(input$module)
                    fig_eigengene(wg, input$module, groups_vec(), mode) },
        network = { shiny::req(input$module)
                    fig_module_network(wg, input$module, mode = mode) })
    })

    output$plot <- shiny::renderPlot({ print(cur_plot()) })

    output$hub_table <- DT::renderDT({
      wg <- wg_rv(); shiny::req(wg, input$module)
      DT::datatable(hub_genes(wg, input$module, n = 25),
                    rownames = FALSE, class = "compact stripe hover",
                    options = list(pageLength = 8, dom = "tp"))
    })

    # Single selected module
    shiny::observeEvent(input$enrich, {
      wg <- wg_rv(); shiny::req(wg, input$module)
      org <- organism_reactive() %||% "human"
      genes <- names(wg$modules)[wg$modules == input$module]
      enrich_view("single")
      shiny::withProgress(message = "Enriching module...", value = 0.4, {
        tryCatch({
          enrich_rv(run_ora(genes, org, db = "GO", ont = "BP",
                            universe = colnames(wg$datExpr)))
          if (nrow(enrich_rv()) == 0)
            shiny::showNotification("No enriched GO BP terms for this module.",
                                    type = "warning", duration = 5)
        }, error = function(e) {
          enrich_rv(NULL)
          shiny::showNotification(paste("Module enrichment failed:",
                                        conditionMessage(e)),
                                  type = "error", duration = 10)
        })
      })
    })

    # All modules at once (compareCluster-style)
    shiny::observeEvent(input$enrich_all, {
      wg <- wg_rv(); shiny::req(wg)
      org <- organism_reactive() %||% "human"
      enrich_view("all")
      shiny::withProgress(message = "Enriching all modules...", value = 0.3, {
        tryCatch({
          modules_rv(enrich_modules(wg, org, db = "GO", ont = "BP", n_per = 5))
        }, error = function(e) {
          modules_rv(NULL)
          shiny::showNotification(paste("Module enrichment failed:",
                                        conditionMessage(e)),
                                  type = "error", duration = 10)
        })
      })
    })

    output$enrich_plot <- shiny::renderPlot({
      if (enrich_view() == "all") {
        shiny::req(modules_rv())
        print(fig_module_enrichment(modules_rv()))
      } else {
        er <- enrich_rv(); shiny::req(er, nrow(er) > 0)
        print(fig_enrich_dot(er, n = 15))
      }
    })

    output$enrich_table <- DT::renderDT({
      df <- if (enrich_view() == "all") {
        shiny::req(modules_rv())
        modules_rv()[, intersect(c("module", "Description", "Count", "padj"),
                                 colnames(modules_rv()))]
      } else {
        er <- enrich_rv(); shiny::req(er, nrow(er) > 0)
        er[, intersect(c("Description", "Count", "padj"), colnames(er))]
      }
      df$padj <- formatC(df$padj, format = "e", digits = 2)
      DT::datatable(df, rownames = FALSE, class = "compact stripe hover",
                    options = list(pageLength = 8, dom = "tp"))
    })

    output$net_dl <- shiny::downloadHandler(
      filename = function() paste0("wgcna_", input$view %||% "plot", ".",
                                   input$net_fmt),
      content = function(file) {
        save_ggplot(cur_plot(), file, input$net_fmt,
                    input$net_w, input$net_h, as.integer(input$net_dpi))
      }
    )
  })
}
