#' Project state management
#'
#' Save and restore complete RNAflow analysis sessions. A project bundles
#' the input data, all parameters, and results into a single .rds file
#' that can be reopened later or shared with collaborators.
#'
#' @name project_state
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
    de_results  = NULL,                 # data.frame from DESeq2
    de_params   = list(),               # design formula, contrast, etc.
    figures     = list(),               # cached plot params, not the plots themselves
    enrichment  = list(),               # GSEA/ORA results
    wgcna       = list(),               # WGCNA modules + correlations
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
  obj
}
