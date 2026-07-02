#' Reusable presentational UI components
#'
#' Small, dependency-light building blocks that give every tab a consistent
#' look: status banners, empty-state placeholders, page headers, and stat
#' tiles. These are pure view helpers -- they return \pkg{htmltools} tags and
#' carry no server logic, so they are safe to drop into any module UI. Styling
#' lives in \code{inst/app/www/rnaflow.css}.
#'
#' @name ui_components
#' @keywords internal
NULL

#' Status banner (info / warning / danger / success)
#'
#' A compact, consistently-styled callout with an optional leading icon.
#'
#' @param ... banner content (text or tags)
#' @param type one of "info", "warning", "danger", "success"
#' @param icon optional Font Awesome icon name; a sensible default is chosen
#'   per type when `NULL`. Pass `FALSE` to omit the icon.
#' @return a \code{div} tag
#' @keywords internal
ui_banner <- function(..., type = c("info", "warning", "danger", "success"),
                      icon = NULL) {
  type <- match.arg(type)
  cls <- switch(type, info = "", warning = " rf-warning",
                danger = " rf-danger", success = " rf-success")
  ic <- NULL
  if (!identical(icon, FALSE)) {
    nm <- icon %||% switch(type,
                           info = "circle-info", warning = "triangle-exclamation",
                           danger = "circle-xmark", success = "circle-check")
    ic <- shiny::icon(nm, class = "rf-ic")
  }
  shiny::div(class = paste0("rnaflow-banner", cls), ic, shiny::span(...))
}

#' Empty-state placeholder
#'
#' Professional guidance to show in a tab before the required inputs exist
#' (no data, no DE results, no enrichment run yet), instead of a raw blank or
#' an error.
#'
#' @param title short headline (e.g. "No differential-expression results yet")
#' @param message one or more guidance sentences (text or tags)
#' @param icon Font Awesome icon name (default "inbox")
#' @return a \code{div} tag
#' @keywords internal
ui_empty_state <- function(title, message = NULL, icon = "inbox") {
  shiny::div(
    class = "rnaflow-empty",
    shiny::icon(icon, class = "rf-empty-ic"),
    shiny::div(class = "rf-empty-title", title),
    if (!is.null(message)) shiny::div(class = "rf-empty-msg", message)
  )
}

#' Tab page header with optional microcopy
#'
#' @param title the tab's title
#' @param subtitle optional one-line description shown under the title
#' @return a \code{div} tag
#' @keywords internal
ui_page_header <- function(title, subtitle = NULL) {
  shiny::div(
    class = "rnaflow-page-header",
    shiny::div(class = "rf-h", title),
    if (!is.null(subtitle)) shiny::div(class = "rf-sub", subtitle)
  )
}

#' Coloured stat tile (big number + label)
#'
#' @param n the number / value to show
#' @param lbl the caption under it
#' @param bg background colour (hex)
#' @return a \code{div} tag
#' @keywords internal
ui_stat_tile <- function(n, lbl, bg = "#1D9E75") {
  shiny::div(
    class = "rnaflow-tile", style = paste0("background:", bg, ";"),
    shiny::div(class = "stat-n", n),
    shiny::div(class = "stat-lbl", lbl)
  )
}
