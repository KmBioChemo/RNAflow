#' Activity figures
#'
#' Diverging bar chart of the top transcription-factor / pathway activity
#' scores from [run_activity()]. Pure ggplot, using the shared RNAflow theme.
#'
#' @name fig_decoupler
#' @keywords internal
NULL

#' Diverging bar chart of top activity scores
#'
#' @param activity an activity data.frame (`source`, `score`, `padj`) from
#'   [run_activity()]
#' @param n number of top regulators / pathways (by |score|) to show
#' @param col_up,col_down colors for activated / repressed bars
#' @param title optional plot title
#' @param mode "exploration" or "publication" (figure theme)
#' @return a ggplot object
#' @export
fig_activity_bar <- function(activity, n = 20, col_up = "#C0392B",
                             col_down = "#2980B9", title = NULL,
                             mode = c("exploration", "publication")) {
  mode <- match.arg(mode)
  if (!is.data.frame(activity) || nrow(activity) == 0) {
    stop("No activity scores to plot.", call. = FALSE)
  }
  d <- activity[order(-abs(activity$score)), , drop = FALSE]
  d <- utils::head(d, n)
  d <- d[order(d$score), , drop = FALSE]
  d$source <- factor(d$source, levels = d$source)
  d$dir <- ifelse(d$score >= 0, "Activated", "Repressed")
  ggplot2::ggplot(d, ggplot2::aes(x = .data$score, y = .data$source,
                                  fill = .data$dir)) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.3, color = "grey55") +
    ggplot2::scale_fill_manual(
      values = c("Activated" = col_up, "Repressed" = col_down), name = NULL) +
    ggplot2::labs(x = "Activity score (t-value)", y = NULL, title = title) +
    fig_theme(mode)
}
