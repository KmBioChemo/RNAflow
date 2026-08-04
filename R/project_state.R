#' Project state management
#'
#' Save and restore complete RNAflow analysis sessions. A project bundles
#' the input data, all parameters, and results into a single .rds file
#' that can be reopened later or shared with collaborators.
#'
#' @name project_state
#' @keywords internal
NULL

#' Create an empty project state
#'
#' @param name project name (used as default filename)
#' @return a list with the canonical project structure
#' @keywords internal
empty_project <- function(name = "untitled") {
  list(
    name        = name,
    created_at  = Sys.time(),
    modified_at = Sys.time(),
    rnaflow_version = utils::packageVersion("RNAflow"),
    organism    = NA_character_,        # "human" / "mouse" / "rat"
    counts      = NULL,                 # numeric matrix
    metadata    = NULL,                 # data.frame
    de_results  = NULL,                 # data.frame from DESeq2 (active contrast)
    de_params   = list(),               # design formula, contrast, etc.
    normalization_method = "vst",       # "vst" or "log2(counts+1) [VST fallback]"
    contrasts   = list(),               # named contrast store (multi-contrast)
    figures     = list(),               # cached plot params, not the plots themselves
    enrichment  = list(),               # GSEA/ORA results
    wgcna       = list(),               # WGCNA modules + correlations
    activity    = list(),               # decoupleR TF/pathway activity (type, table)
    signatures  = list(),               # GSVA/ssGSEA per-sample scores (settings + matrix)
    ai_interpretation = NULL,           # optional AI narrative (list: text, model)
    notes       = character(0)          # free-text annotations
  )
}

#' Save a project to disk
#'
#' @param project a project list (from [empty_project()] or session)
#' @param path file path (will get .rnaflow.rds extension if missing)
#' @return invisibly returns the path written to
#' @export
save_project <- function(project, path) {
  if (!grepl("\\.rnaflow\\.rds$", path)) {
    path <- paste0(sub("\\.rds$", "", path), ".rnaflow.rds")
  }
  project$modified_at <- Sys.time()
  saveRDS(project, path)
  invisible(path)
}

#' Load a project from disk
#'
#' @param path file path
#' @return the project list
#' @export
load_project <- function(path) {
  if (!file.exists(path)) {
    stop("Project file not found: ", path, call. = FALSE)
  }
  obj <- tryCatch(readRDS(path),
                  error = function(e) stop("Could not read project: ",
                                           conditionMessage(e), call. = FALSE))
  if (!is.list(obj) || !"rnaflow_version" %in% names(obj)) {
    stop("File does not look like a valid RNAflow project.", call. = FALSE)
  }
  # Backfill any slots added in newer versions so projects saved by older
  # RNAflow releases load with the canonical structure (missing fields take
  # their empty-project defaults rather than being absent).
  defaults <- empty_project(obj$name %||% "untitled")
  for (k in setdiff(names(defaults), names(obj))) obj[[k]] <- defaults[[k]]
  obj
}

#' Assemble a project from the current session state
#'
#' Bundles the live analysis objects into the canonical project structure
#' (see [empty_project()]). Shared by the project-manager and report modules.
#'
#' @param name project name
#' @param organism organism keyword
#' @param counts counts matrix (or NULL)
#' @param metadata metadata data.frame (or NULL)
#' @param contrasts the contrast store (named list)
#' @param settings optional list with `enrichment` / `wgcna` parameter records
#'   (captured by the Enrichment / Network tabs) for exact reproducibility
#' @return a project list
#' @export
assemble_project <- function(name, organism = NA_character_,
                             counts = NULL, metadata = NULL,
                             contrasts = list(), settings = list()) {
  p <- empty_project(if (!is.null(name) && nzchar(name)) name else "untitled")
  p$organism  <- organism %||% NA_character_
  p$counts    <- counts
  p$metadata  <- metadata
  p$contrasts <- contrasts %||% list()
  p$enrichment <- if (is.null(settings$enrichment)) list() else settings$enrichment
  p$wgcna      <- if (is.null(settings$wgcna)) list() else settings$wgcna
  p$activity   <- if (is.null(settings$activity)) list() else settings$activity
  # Signatures (GSVA/ssGSEA): canonical slot is `signatures`; fall back to the
  # legacy `gsva` key for sessions/records created before 0.14.1.
  p$signatures <- settings$signatures %||% settings$gsva %||% list()
  p$ai_interpretation <- settings$ai_interpretation
  p$normalization_method <- settings$normalization_method %||% "vst"
  active <- contrast_store_results(p$contrasts)
  p$de_results <- if (length(active)) active[[1]] else NULL
  p
}

#' Insert or update a contrast in a contrast store
#'
#' A contrast store is a named list keyed by contrast label. Each entry holds
#' the DE results data.frame, the parameters used to compute it, and a
#' timestamp. Re-adding an existing label updates that entry in place
#' (keeping its position), so re-running the same contrast refreshes it
#' rather than duplicating it.
#'
#' @param store the current store (a named list; may be empty)
#' @param label contrast label (unique key)
#' @param results DE results data.frame
#' @param params named list of parameters used to compute the contrast
#' @param created optional timestamp (defaults to now)
#' @return the updated store
#' @export
contrast_store_upsert <- function(store, label, results, params = list(),
                                  created = Sys.time()) {
  if (is.null(store)) store <- list()
  label <- as.character(label)
  if (!nzchar(label)) stop("Contrast label must be non-empty.", call. = FALSE)
  store[[label]] <- list(
    label   = label,
    results = results,
    params  = params,
    created = created
  )
  store
}

#' Directory where recent projects are cached
#'
#' Defaults to a per-user data directory (`tools::R_user_dir`), overridable
#' via `options(rnaflow.recent_dir = ...)` (used in tests). Created on first
#' use.
#'
#' @return the directory path (created if missing)
#' @keywords internal
rnaflow_recent_dir <- function() {
  dir <- getOption("rnaflow.recent_dir",
                   tools::R_user_dir("RNAflow", "data"))
  dir <- file.path(dir, "projects")
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  dir
}

# Small stable non-cryptographic hash (base R only) so distinct project names
# that sanitise to the same string don't collide in the recent cache, while the
# same name still maps to one file (idempotent re-caching).
name_hash <- function(s) {
  h <- 0
  for (ch in utf8ToInt(enc2utf8(s))) h <- (h * 33 + ch) %% 2147483647
  sprintf("%09.0f", h)
}

#' Cache a project file in the recent-projects directory
#'
#' Copies a saved/loaded `.rnaflow.rds` into the recent directory so it can
#' be re-opened from the launch panel later.
#'
#' @param path path to an existing project file
#' @param name optional display name used to build the cached filename
#' @return the cached file path (invisibly)
#' @keywords internal
cache_recent_project <- function(path, name = NULL) {
  if (!file.exists(path)) return(invisible(NULL))
  base <- if (!is.null(name) && nzchar(name)) {
    paste0(gsub("[^A-Za-z0-9_-]+", "_", name), "-", name_hash(name))
  } else {
    sub("\\.rnaflow\\.rds$", "", basename(path))
  }
  dest <- file.path(rnaflow_recent_dir(), paste0(base, ".rnaflow.rds"))
  file.copy(path, dest, overwrite = TRUE)
  invisible(dest)
}

#' List recent cached projects
#'
#' @param max_n maximum number of entries to return (most recent first)
#' @return a data.frame with columns `file`, `name`, `modified_at` (possibly
#'   zero rows)
#' @keywords internal
list_recent_projects <- function(max_n = 8) {
  dir <- rnaflow_recent_dir()
  files <- list.files(dir, pattern = "\\.rnaflow\\.rds$", full.names = TRUE)
  if (length(files) == 0) {
    return(data.frame(file = character(0), name = character(0),
                      modified_at = as.POSIXct(character(0))))
  }
  rows <- lapply(files, function(f) {
    info <- tryCatch(readRDS(f), error = function(e) NULL)
    nm <- if (is.list(info) && !is.null(info$name)) info$name
          else sub("\\.rnaflow\\.rds$", "", basename(f))
    mt <- if (is.list(info) && !is.null(info$modified_at)) info$modified_at
          else file.info(f)$mtime
    data.frame(file = f, name = nm, modified_at = as.POSIXct(mt),
               stringsAsFactors = FALSE)
  })
  df <- do.call(rbind, rows)
  df <- df[order(df$modified_at, decreasing = TRUE), , drop = FALSE]
  utils::head(df, max_n)
}

#' Extract the bare DE results from a contrast store
#'
#' @param store a contrast store (named list)
#' @return a named list of DE results data.frames (the shape expected by the
#'   [analysis_compare] and [fig_compare] functions)
#' @keywords internal
contrast_store_results <- function(store) {
  if (is.null(store) || length(store) == 0) return(list())
  stats::setNames(lapply(store, function(e) e$results), names(store))
}
