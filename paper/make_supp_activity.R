# Supplementary figure S1 -- pathway activity inference (PROGENy) on the airway
# dexamethasone contrast. This is the *demonstrated* half of the activity module:
# PROGENy recovers the expected glucocorticoid / anti-inflammatory pathway signal.
# (TF activity via CollecTRI stays available in the app but is not featured in the
# manuscript -- its live OmniPath fetch is not reliable enough to guarantee for a
# reviewer, and this figure must be reproducible offline.)
#
# NOTE on the call: PROGENy networks carry a `weight` column and must be scored
# with the multivariate model, i.e. run_activity(..., method = "mlm",
# mor_col = "weight"). The default (method = "ulm", mor_col = "mor") is for
# CollecTRI TF regulons and errors on a PROGENy network -- that mismatch is why
# the old gallery panel `18_pathway_activity` never rendered.
#
# The PROGENy network is fetched from OmniPath ONCE and cached to
# paper/.activity_cache.rds, so every subsequent run (and CI, and a reviewer's
# re-run) is fully offline and reproducible. Delete the cache to refetch.
#
# Run from the package root, on a machine where OmniPath is reachable the first
# time:  Rscript paper/make_supp_activity.R
suppressPackageStartupMessages({
  if (!suppressWarnings(require(RNAflow, quietly = TRUE))) devtools::load_all(".", quiet = TRUE)
  library(ggplot2)
})
set.seed(1)
ext   <- function(f) system.file("extdata", f, package = "RNAflow")
cache <- "paper/.activity_cache.rds"

## ---- compute once, cache (network fetch is the only online step) -------
if (file.exists(cache)) {
  D <- readRDS(cache); message("loaded activity cache")
} else {
  message("computing (fetches the PROGENy network from OmniPath once)...")
  ac <- read_counts(ext("demo_airway_counts.csv"))
  am <- read_metadata(ext("demo_airway_metadata.csv"), counts_samples = colnames(ac))
  a_res <- run_deseq2(ac, am, design = ~ cell + condition,
                      contrast = c("condition", "Dex", "Control"), shrink = TRUE)
  # PROGENy pathway network (weight column) -> multivariate linear model.
  progeny_net <- get_pathway_network("human", top = 500)
  act <- run_activity(a_res, progeny_net, method = "mlm",
                      mor_col = "weight", by = "stat", min_size = 5)
  D <- list(act = act, net_rows = nrow(progeny_net))
  saveRDS(D, cache)
  message(sprintf("cached: %d pathways scored from a %d-edge PROGENy network",
                  nrow(D$act), D$net_rows))
}

## ---- render the supplementary figure ----------------------------------
p <- fig_activity_bar(D$act, n = 14, mode = "publication") +
  labs(title = "Pathway activity (PROGENy) — dexamethasone vs control (airway)") +
  theme(plot.title = element_text(size = 11, face = "bold", margin = margin(b = 4)))

ggsave("paper/figures/figureS1_activity.png", p, width = 7.5, height = 5,
       dpi = 300, bg = "white")
ggsave("paper/figures/figureS1_activity.pdf", p, width = 7.5, height = 5, bg = "white")
cat("wrote figureS1_activity\n")

## ---- report the top pathways so the caption/text can quote real numbers
top <- head(D$act[order(-abs(D$act$score)), ], 6)
cat("\nTop pathways by |activity score| (fill these into paper.md TODO):\n")
print(top[, intersect(c("source", "score", "p_value", "padj"), colnames(top))],
      row.names = FALSE)
