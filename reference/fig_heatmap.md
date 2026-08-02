# Heatmap figures

Differential expression heatmap with z-scored rows, sample annotations,
and optional direction (up/down) annotation. Requires a real counts
matrix – the silent simulation from the original app has been removed to
prevent accidentally publishing fake data.

## Usage

``` r
fig_heatmap(
  counts_mat,
  res,
  metadata,
  group_colors = list(),
  n_genes = 40,
  gene_src = c("top_n", "all_sig"),
  padj_thr = 0.05,
  lfc_thr = 1,
  palette_name = "RdBu",
  title = NULL,
  show_title = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  show_dend_rows = TRUE,
  show_dend_cols = TRUE,
  show_legend = TRUE,
  annotation_legend = TRUE,
  show_annotation_names = TRUE,
  direction_annotation = FALSE,
  col_dir_up = "#C0392B",
  col_dir_down = "#2980B9",
  ann_title = ""
)
```

## Arguments

- counts_mat:

  normalized or raw counts matrix (used for expression values)

- res:

  DE results data.frame (used to pick genes by significance)

- metadata:

  sample metadata (column 1 = sample, column 2 = annotation)

- group_colors:

  named list of user-supplied colors keyed by
  `paste0("grp_col_", sanitized_condition)`

- n_genes:

  number of top genes to display (when gene_src == "top_n")

- gene_src:

  "top_n" (default) or "all_sig" (all genes passing thresholds)

- padj_thr, lfc_thr:

  significance thresholds

- palette_name:

  colorbrewer / viridis palette name

- title:

  character; title of the heatmap

- show_title:

  whether to draw the title

- cluster_rows, cluster_cols:

  cluster genes / samples

- show_rownames, show_colnames:

  display row / column labels

- show_dend_rows, show_dend_cols:

  display dendrograms

- show_legend, annotation_legend:

  display legends

- show_annotation_names:

  show the row annotation track name (e.g. "Direction") beside the
  track; when FALSE it appears only in the legend. The column annotation
  name is never drawn on the plot (it duplicates the annotation legend
  title and can overflow onto the legends).

- direction_annotation:

  add up/down row annotation

- col_dir_up, col_dir_down:

  colors for direction annotation

- ann_title:

  custom column annotation header (e.g. "Group")

## Value

a pheatmap object (gtable)
