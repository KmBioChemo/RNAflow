# Diverging bar chart of top activity scores

Diverging bar chart of top activity scores

## Usage

``` r
fig_activity_bar(
  activity,
  n = 20,
  col_up = "#C0392B",
  col_down = "#2980B9",
  title = NULL,
  mode = c("exploration", "publication")
)
```

## Arguments

- activity:

  an activity data.frame (`source`, `score`, `padj`) from
  [`run_activity()`](https://KmBioChemo.github.io/RNAflow/reference/run_activity.md)

- n:

  number of top regulators / pathways (by \|score\|) to show

- col_up, col_down:

  colors for activated / repressed bars

- title:

  optional plot title

- mode:

  "exploration" or "publication" (figure theme)

## Value

a ggplot object
