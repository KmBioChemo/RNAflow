# Driver for the network-restricted env: local Hallmark sets + enricher-based ORA,
# so make_panels.R / make_validation.R render without CRAN/Bioc annotation downloads.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
devtools::load_all(".", quiet = TRUE)
read_gmt <- function(f){ l <- strsplit(readLines(f), "\t")
  setNames(lapply(l, function(x) x[-(1:2)]), vapply(l, `[`, "", 1)) }
SETS <- read_gmt("paper/genesets/h.all.v2023.2.Hs.symbols.gmt")
T2G  <- do.call(rbind, lapply(names(SETS), function(n) data.frame(term = n, gene = SETS[[n]])))
# global overrides shadow the package versions when called from the script
get_gene_sets <- function(...) SETS
# run_ora and enrich_modules now use the package's NATIVE GO implementation
# (real GO.db + org.Hs.eg.db are installed) -> real GO-BP panels.
cat(">>> overrides in place; sourcing make_panels.R\n")
source("paper/make_panels.R")
cat(">>> sourcing make_validation.R\n")
tryCatch(source("paper/make_validation.R"), error = function(e) cat("VALIDATION_ERR:", conditionMessage(e), "\n"))
cat(">>> ALL_DONE\n")
