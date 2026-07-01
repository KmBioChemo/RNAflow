# Cross-contrast log2FoldChange signature heatmap

Genes x contrasts heatmap of log2FoldChange, for spotting shared vs.
contrast-specific transcriptional signatures. Genes absent from a
contrast are filled with 0 (no change) so clustering stays well-defined.

## Usage

``` r
fig_lfc_heatmap(
  contrasts,
  genes = NULL,
  gene_src = c("sig_union", "top_var"),
  n_genes = 50,
  padj_thr = 0.05,
  lfc_thr = 1,
  palette_name = "RdBu",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = NULL,
  title = NULL,
  show_title = TRUE
)
```

## Arguments

- contrasts:

  a named list of DE results data.frames

- genes:

  optional character vector to display (overrides gene_src)

- gene_src:

  "sig_union" (default; union of significant genes across contrasts) or
  "top_var" (most variable genes among all)

- n_genes:

  max number of genes to display (top by cross-contrast variance when
  more are available)

- padj_thr, lfc_thr:

  thresholds used when `gene_src == "sig_union"`

- palette_name:

  diverging palette name (default "RdBu")

- cluster_rows, cluster_cols:

  cluster genes / contrasts

- show_rownames:

  display gene labels (default: auto, when \<= 60 genes)

- title, show_title:

  heatmap title

## Value

a pheatmap object (gtable)
