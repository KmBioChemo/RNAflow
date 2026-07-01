# Run DESeq2 on a counts matrix

Constructs a DESeqDataSet from counts and metadata, runs DESeq(), and
extracts results for a user-specified contrast. Optionally applies LFC
shrinkage (apeglm by default).

## Usage

``` r
run_deseq2(
  counts,
  metadata,
  design = ~condition,
  contrast = NULL,
  shrink = TRUE,
  shrink_type = c("apeglm", "ashr", "normal"),
  min_count = 10,
  alpha = 0.05
)
```

## Arguments

- counts:

  validated counts matrix (genes x samples, integer)

- metadata:

  data.frame with sample ID in column 1

- design:

  a one-sided formula referring to columns of `metadata`, e.g.
  `~ condition` or `~ batch + condition`. The variable of interest
  should be the last term.

- contrast:

  a length-3 character vector: c(variable, level_treated,
  level_reference). Example: `c("condition", "Treatment", "Control")`.

- shrink:

  logical, apply LFC shrinkage to the effect-size estimate (recommended
  for ranking / visualization). Inference is unaffected.

- shrink_type:

  shrinkage estimator: "apeglm" (default), "ashr", or "normal"

- min_count:

  minimum row sum to keep a gene (default 10)

- alpha:

  FDR threshold passed to `results()` for independent filtering

## Value

a tidy data.frame with columns: gene, baseMean, log2FoldChange, lfcSE,
stat, pvalue, padj. The estimator actually used is recorded in
`attr(result, "shrink")`.

## Details

Shrinkage adjusts only the effect-size estimates (`log2FoldChange`,
`lfcSE`), improving ranking and visualization. Inference is unaffected:
`stat`, `pvalue` and `padj` always come from the unshrunken Wald test.
Consequently the default GSEA ranking (`rank_genes(by = "stat")`) uses
the unshrunken Wald statistic even when `shrink = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
  counts <- read_counts("counts.csv")
  meta   <- read_metadata("metadata.csv", counts_samples = colnames(counts))
  res    <- run_deseq2(counts, meta,
                       design = ~ condition,
                       contrast = c("condition", "Treatment", "Control"))
} # }
```
