#' Reusable UI widgets
#'
#' Shared components used across RNAflow modules: color pickers with
#' synchronized text fields, slider+numeric input pairs, export toolbar.
#'
#' @name ui_widgets
#' @keywords internal
NULL

#' Section title (small caps style)
#' @param lbl text to display
#' @keywords internal
ui_section_title <- function(lbl) {
  shiny::tags$span(lbl, class = "s-title")
}

#' Color picker with synchronized hex text input
#'
#' @param id input ID (Shiny will see input[[id]] with the hex value)
#' @param lbl display label
#' @param default starting hex color
#' @return a UI element
#' @keywords internal
ui_color_picker <- function(id, lbl, default) {
  js_cp <- sprintf(
    "Shiny.setInputValue('%s', this.value, {priority:'event'});\
    document.getElementById('txt_%s').value=this.value;", id, id
  )
  js_tx <- sprintf(
    "if(/^#[0-9A-Fa-f]{6}$/.test(this.value)){\
    Shiny.setInputValue('%s', this.value, {priority:'event'});\
    document.getElementById('cp_%s').value=this.value;}", id, id
  )
  shiny::div(
    style = "display:flex;align-items:center;gap:4px;margin-bottom:3px;",
    shiny::tags$label(lbl, style = "min-width:60px;margin:0;white-space:nowrap;"),
    shiny::tags$input(
      id = paste0("cp_", id), type = "color", value = default,
      style = "width:24px;height:22px;padding:1px;border:1px solid #ced4da;border-radius:3px;cursor:pointer;flex-shrink:0;",
      oninput = js_cp
    ),
    shiny::tags$input(
      id = paste0("txt_", id), type = "text", value = default,
      placeholder = "#RRGGBB", class = "form-control",
      style = "font-family:monospace;height:22px;padding:2px 5px;",
      oninput = js_tx
    )
  )
}

#' Synchronized slider + numeric input pair
#'
#' Useful for fine control when a slider alone is too coarse.
#'
#' @param slider_id slider input ID
#' @param num_id numeric input ID
#' @param lbl label
#' @param mn,mx min/max
#' @param val initial value
#' @param stp step
#' @return a UI element
#' @keywords internal
ui_slider_num <- function(slider_id, num_id, lbl, mn, mx, val, stp) {
  shiny::tagList(
    shiny::tags$label(lbl),
    shiny::div(
      style = "display:flex;align-items:center;gap:3px;",
      shiny::div(style = "flex:1;",
                 shiny::sliderInput(slider_id, NULL, min = mn, max = mx,
                                    value = val, step = stp, ticks = FALSE)),
      shiny::div(style = "width:62px;",
                 shiny::numericInput(num_id, NULL, val, min = mn, max = mx,
                                     step = stp, width = "100%"))
    )
  )
}

#' Export toolbar (format + dimensions + download button)
#'
#' @param prefix prefix used to build input IDs (`<prefix>_fmt`, `_w`, `_h`,
#'   `_dpi`, `_dl`)
#' @param def_w,def_h default dimensions (inches)
#' @return a UI element
#' @keywords internal
ui_export_bar <- function(prefix, def_w, def_h) {
  shiny::div(
    style = "border-top:1px solid #e9ecef;margin-top:5px;padding-top:5px;",
    ui_section_title("Export"),
    shiny::div(
      style = "display:flex;gap:3px;flex-wrap:wrap;align-items:flex-end;",
      shiny::div(
        style = "width:70px;", shiny::tags$label("Format"),
        shiny::selectInput(paste0(prefix, "_fmt"), NULL,
                           choices = c("PDF" = "pdf", "PNG" = "png", "TIFF" = "tiff"),
                           selected = "pdf", width = "100%")
      ),
      shiny::div(style = "width:52px;", shiny::tags$label("W (in)"),
                 shiny::numericInput(paste0(prefix, "_w"), NULL, def_w,
                                     min = 2, max = 30, step = 0.5, width = "100%")),
      shiny::div(style = "width:52px;", shiny::tags$label("H (in)"),
                 shiny::numericInput(paste0(prefix, "_h"), NULL, def_h,
                                     min = 2, max = 30, step = 0.5, width = "100%")),
      shiny::div(style = "width:54px;", shiny::tags$label("DPI"),
                 shiny::selectInput(paste0(prefix, "_dpi"), NULL,
                                    choices = c("150", "300", "600"),
                                    selected = "300", width = "100%")),
      shiny::downloadButton(paste0(prefix, "_dl"), "Export",
                            class = "btn btn-sm btn-success mt-1")
    )
  )
}

#' Stat badge (number + label, colored background)
#' @keywords internal
ui_stat_badge <- function(n, lbl, bg) {
  shiny::div(
    class = "stat-badge", style = paste0("background:", bg, ";"),
    shiny::div(class = "stat-n", n),
    shiny::div(class = "stat-lbl", lbl)
  )
}

#' Advanced settings accordion wrapper
#' @keywords internal
ui_advanced_panel <- function(...) {
  bslib::accordion(
    class = "adv-acc", open = FALSE,
    bslib::accordion_panel("Advanced settings", ...)
  )
}
