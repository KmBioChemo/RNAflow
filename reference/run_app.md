# Launch the RNAflow Shiny application

Runs the full RNAflow Shiny app, which assembles all modules (data, DE,
volcano, heatmap, PCA, ...) into a single interface.

## Usage

``` r
run_app(port = NULL, launch_browser = TRUE)
```

## Arguments

- port:

  port to launch on (default: random free port)

- launch_browser:

  if TRUE, open in default browser

## Value

invisibly returns the Shiny app object

## Examples

``` r
if (FALSE) { # \dontrun{
  RNAflow::run_app()
} # }
```
