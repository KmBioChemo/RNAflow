#' Package hooks
#'
#' Registers the bundled web assets (`inst/app/www`) under a stable resource
#' prefix so the Shiny UI can link to the stylesheet with `href = "rnaflow/..."`.
#' Without this the `<link>` in [app_ui()] would 404 and none of the app's
#' styling would load.
#'
#' @name rnaflow-package-hooks
#' @keywords internal
NULL

.onLoad <- function(libname, pkgname) {
  www <- system.file("app/www", package = pkgname)
  if (nzchar(www) && dir.exists(www)) {
    shiny::addResourcePath("rnaflow", www)
  }
  invisible()
}

.onUnload <- function(libpath) {
  suppressWarnings(try(shiny::removeResourcePath("rnaflow"), silent = TRUE))
  invisible()
}
