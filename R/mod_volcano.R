#' Volcano plot module
#'
#' Shiny module wrapping [fig_volcano()] and [fig_volcano_interactive()],
#' with all the controls from the original app.
#'
#' @param id namespace ID
#' @name mod_volcano
NULL

#' @rdname mod_volcano
#' @export
mod_volcano_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::layout_sidebar(
    sidebar = bslib::sidebar(
      width = 320,
      ui_section_title("Thresholds"),
      ui_slider_num(ns("lfc_sld"), ns("lfc_num"),
                    "|log2 fold change| >", 0, 5, 1, 0.1),
      ui_slider_num(ns("padj_sld"), ns("padj_num"),
                    "FDR (padj) <", 0.001, 0.2, 0.05, 0.001),
      ui_section_title("Appearance"),
      ui_slider_num(ns("pts_sld"), ns("pts_num"),
                    "Point size", 0.5, 5, 1.8, 0.1),
      ui_slider_num(ns("nlab_sld"), ns("nlab_num"),
                    "Top genes to label", 0, 80, 20, 1),
      shiny::selectInput(ns("leg"), "Legend position",
                         choices = c("Right", "Top-left", "Top-right",
                                     "Bottom-left", "Bottom-right", "None")),
      shiny::checkboxInput(ns("publication"),
                           "Publication mode (export)", value = FALSE),
      ui_advanced_panel(
        ui_color_picker(ns("cu"), "Up",     "#C0392B"),
        ui_color_picker(ns("cd"), "Down",   "#2980B9"),
        ui_color_picker(ns("cn"), "NS",     "#BDC3C7"),
        ui_color_picker(ns("cc"), "Cutoff", "#7F8C8D"),
        shiny::textInput(ns("xlab"), "X-axis label", ""),
        shiny::textInput(ns("ylab"), "Y-axis label", ""),
        shiny::textInput(ns("title"), "Title", ""),
        shiny::checkboxInput(ns("show_title"), "Show title", TRUE),
        shiny::checkboxInput(ns("show_subtitle"), "Show subtitle (counts)", TRUE),
        shiny::numericInput(ns("x_min"), "X min (blank = auto)", NA),
        shiny::numericInput(ns("x_max"), "X max (blank = auto)", NA),
        shiny::numericInput(ns("y_max"), "Y max (blank = auto)", NA)
      ),
      ui_export_bar(ns("vol"), 8, 6)
    ),
    shiny::uiOutput(ns("stats")),
    shinycssloaders::withSpinner(
      plotly::plotlyOutput(ns("plot_interactive"), height = "560px"),
      type = 6, color = "#1D9E75"
    ),
    shiny::tags$hr(),
    shiny::h5("Significant genes"),
    DT::DTOutput(ns("de_table"))
  )
}

#' @rdname mod_volcano
#' @param de_reactive a reactive returning a DE results data.frame
#' @export
mod_volcano_server <- function(id, de_reactive) {
  shiny::moduleServer(id, function(input, output, session) {

    # Mirror slider <-> numeric
    mirror <- function(s, n) {
      shiny::observeEvent(input[[s]],
        shiny::updateNumericInput(session, n, value = input[[s]]),
        ignoreInit = TRUE)
      shiny::observeEvent(input[[n]],
        shiny::updateSliderInput(session, s, value = input[[n]]),
        ignoreInit = TRUE)
    }
    mirror("lfc_sld", "lfc_num")
    mirror("padj_sld", "padj_num")
    mirror("pts_sld", "pts_num")
    mirror("nlab_sld", "nlab_num")

    cu <- shiny::reactive(safe_col(input$cu, "#C0392B"))
    cd <- shiny::reactive(safe_col(input$cd, "#2980B9"))
    cn <- shiny::reactive(safe_col(input$cn, "#BDC3C7"))
    cc <- shiny::reactive(safe_col(input$cc, "#7F8C8D"))

    output$stats <- shiny::renderUI({
      res <- de_reactive()
      shiny::req(res, input$lfc_num, input$padj_num)
      df <- res
      ok <- !is.na(df$padj) & !is.na(df$log2FoldChange)
      df$r <- "NS"
      df$r[ok & df$padj < input$padj_num & df$log2FoldChange >  input$lfc_num] <- "Up"
      df$r[ok & df$padj < input$padj_num & df$log2FoldChange < -input$lfc_num] <- "Down"
      shiny::div(style = "display:flex;gap:7px;margin-bottom:7px;flex-wrap:wrap;",
        ui_stat_badge(sum(df$r == "Up"),   "Up-regulated",    cu()),
        ui_stat_badge(sum(df$r == "Down"), "Down-regulated",  cd()),
        ui_stat_badge(sum(df$r == "NS"),   "Not significant", "#95A5A6"),
        ui_stat_badge(nrow(df),            "Total genes",     "#2C3E50")
      )
    })

    output$plot_interactive <- plotly::renderPlotly({
      res <- de_reactive()
      shiny::req(res, input$lfc_num, input$padj_num)
      fig_volcano_interactive(
        res,
        lfc_thr = input$lfc_num, padj_thr = input$padj_num,
        col_up = cu(), col_down = cd(), col_ns = cn(), col_cut = cc(),
        pt_size = input$pts_num,
        xlab = input$xlab, ylab = input$ylab, title = input$title,
        show_title = input$show_title, show_subtitle = input$show_subtitle,
        leg_pos = input$leg,
        x_min = input$x_min, x_max = input$x_max, y_max = input$y_max
      )
    })

    v_static <- shiny::reactive({
      res <- de_reactive()
      shiny::req(res, input$lfc_num, input$padj_num)
      fig_volcano(
        res,
        lfc_thr = input$lfc_num, padj_thr = input$padj_num,
        n_label = max(0L, as.integer(input$nlab_num)),
        col_up = cu(), col_down = cd(), col_ns = cn(), col_cut = cc(),
        pt_size = input$pts_num,
        xlab = input$xlab, ylab = input$ylab, title = input$title,
        show_title = input$show_title, show_subtitle = input$show_subtitle,
        leg_pos = input$leg,
        x_min = input$x_min, x_max = input$x_max, y_max = input$y_max,
        mode = if (isTRUE(input$publication)) "publication" else "exploration"
      )
    })

    output$vol_dl <- shiny::downloadHandler(
      filename = function() paste0("volcano_plot.", input$vol_fmt),
      content  = function(file) {
        save_ggplot(v_static(), file, input$vol_fmt,
                    input$vol_w, input$vol_h, as.integer(input$vol_dpi))
      }
    )

    output$de_table <- DT::renderDT({
      res <- de_reactive()
      shiny::req(res, input$lfc_num, input$padj_num)
      df <- res
      ok <- !is.na(df$padj) & !is.na(df$log2FoldChange)
      df$regulation <- "NS"
      df$regulation[ok & df$padj < input$padj_num & df$log2FoldChange >  input$lfc_num] <- "Up"
      df$regulation[ok & df$padj < input$padj_num & df$log2FoldChange < -input$lfc_num] <- "Down"
      df <- df[df$regulation != "NS", ]
      df <- df[order(df$padj), ]
      df$log2FoldChange <- round(df$log2FoldChange, 3)
      df$padj <- formatC(df$padj, format = "e", digits = 2)
      DT::datatable(
        df,
        options = list(pageLength = 12, scrollX = TRUE,
                       dom = "<'row'<'col-sm-12't>><'row'<'col-sm-5'i><'col-sm-7'p>>",
                       columnDefs = list(list(className = "dt-center", targets = "_all"))),
        rownames = FALSE, class = "compact stripe hover"
      ) %>%
        DT::formatStyle("regulation", fontWeight = "bold",
                        backgroundColor = DT::styleEqual(
                          c("Up", "Down"), c("#FADBD8", "#D6EAF8")))
    })
  })
}
