#' Color utilities
#'
#' Bulletproof color handling for figure generation. Validates hex codes,
#' falls back gracefully on invalid inputs, and builds n-color vectors that
#' never crash downstream rendering (col2rgb / scale_color_manual).
#'
#' @name color_utils
#' @keywords internal
NULL

#' Hardcoded fallback palette
#'
#' 8 distinct colors guaranteed valid in any context. Used when
#' RColorBrewer or user input fails.
#'
#' @keywords internal
SAFE8 <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
           "#FF7F00", "#A65628", "#F781BF", "#666666")

#' Palette choices for UI dropdowns
#' @keywords internal
PALETTE_CHOICES <- list(
  "Diverging"  = c("RdBu", "RdYlBu", "RdYlGn", "PuOr", "PRGn",
                   "BrBG", "Spectral", "PiYG", "RdGy"),
  "Sequential" = c("Blues", "Reds", "Greens", "Purples", "Oranges",
                   "YlOrRd", "YlGnBu", "BuPu", "GnBu"),
  "Multi-hue"  = c("viridis", "plasma", "magma", "inferno", "cividis")
)

#' Viridis-family palette stops (no dependency on viridisLite)
#' @keywords internal
VIRIDIS_STOPS <- list(
  viridis = c("#440154", "#31688E", "#35B779", "#FDE725"),
  plasma  = c("#0D0887", "#7E03A8", "#CC4678", "#F89441", "#F0F921"),
  magma   = c("#000004", "#3B0F70", "#8C2981", "#DE4968", "#FE9F6D", "#FCFDBF"),
  inferno = c("#000004", "#420A68", "#932667", "#DD513A", "#FCA50A", "#FCFFA4"),
  cividis = c("#00204D", "#31446B", "#666870", "#9F9B63", "#DCD54A", "#FFEA46")
)

#' Test whether a value is a valid 6-digit hex color
#'
#' @param x value to test (typically a string)
#' @return logical of length 1
#' @keywords internal
is_hex6 <- function(x) {
  length(x) == 1 && !is.null(x) && !is.na(x) &&
    grepl("^#[0-9A-Fa-f]{6}$", trimws(as.character(x)))
}

#' Safely coerce a value to a hex color, with fallback
#'
#' Used everywhere we read a user input that should be a color. Guarantees
#' the return value is always a valid hex code, never NA or invalid.
#'
#' @param val candidate value
#' @param fb fallback color if val is invalid
#' @return a valid hex color string
#' @keywords internal
safe_col <- function(val, fb = "#888888") {
  v <- tryCatch(trimws(as.character(val[[1]])), error = function(e) "")
  if (is_hex6(v)) return(v)
  f <- tryCatch(trimws(as.character(fb[[1]])), error = function(e) "#888888")
  if (is_hex6(f)) return(f)
  "#888888"
}

#' Build an n-color palette
#'
#' @param name palette name (from PALETTE_CHOICES or viridis family)
#' @param n number of colors to generate
#' @return character vector of n hex colors
#' @keywords internal
make_palette <- function(name, n = 100) {
  if (name %in% names(VIRIDIS_STOPS)) {
    return(grDevices::colorRampPalette(VIRIDIS_STOPS[[name]])(n))
  }
  raw <- tryCatch(
    RColorBrewer::brewer.pal(RColorBrewer::brewer.pal.info[name, "maxcolors"], name),
    error = function(e) SAFE8
  )
  divs <- c("RdBu", "RdYlBu", "RdYlGn", "PuOr", "PRGn",
            "BrBG", "Spectral", "PiYG", "RdGy")
  grDevices::colorRampPalette(if (name %in% divs) rev(raw) else raw)(n)
}

#' Build n group colors with per-condition override support
#'
#' For each of n conditions, looks up a custom color in `gcl` (a named list
#' from shiny inputs); falls back to a Set1-derived palette if missing or
#' invalid. Guaranteed never to throw an RGB / color parsing error.
#'
#' @param n number of conditions
#' @param conds character vector of condition names (length n)
#' @param gcl named list of user-supplied colors (typically reactive values)
#' @return character vector of n valid hex colors
#' @keywords internal
build_cols <- function(n, conds, gcl) {
  base <- tryCatch(
    grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(max(n, 9)),
    error = function(e) rep(SAFE8, ceiling(n / 8))[seq_len(n)]
  )
  vapply(seq_len(n), function(i) {
    key <- paste0("grp_col_", gsub("[^A-Za-z0-9_]", "_", conds[i]))
    raw <- gcl[[key]]
    col <- safe_col(raw, base[min(i, length(base))])
    tryCatch({ grDevices::col2rgb(col); col },
             error = function(e) SAFE8[((i - 1) %% 8) + 1])
  }, character(1))
}
