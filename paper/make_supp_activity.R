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
  # Prefer the tool's OmniPath fetch; fall back to the offline `progeny` model
  # (same weights, bundled as package data) when OmnipathR/OmniPath is absent.
  progeny_net <- tryCatch(get_pathway_network("human", top = 500), error = function(e) {
    message("OmniPath unavailable (", conditionMessage(e),
            "); building PROGENy network from the offline progeny model.")
    mdl <- progeny::getModel("Human", top = 500)
    do.call(rbind, lapply(colnames(mdl), function(pw) {
      w <- mdl[, pw]; keep <- w != 0
      data.frame(source = pw, target = rownames(mdl)[keep],
                 weight = as.numeric(w[keep]), stringsAsFactors = FALSE)
    }))
  })
  act <- run_activity(a_res, progeny_net, method = "mlm",
                      mor_col = "weight", by = "stat", min_size = 5)
  D <- list(act = act, net_rows = nrow(progeny_net))
  saveRDS(D, cache)
  message(sprintf("cached: %d pathways scored from a %d-edge PROGENy network",
                  nrow(D$act), D$net_rows))
}

## ---- bare panel for the Python plate system (paper/plate/) -----------------
# The composed supplementary figure (figureS1) is assembled by compose.py from
# this bare panel; no standalone plate is written here.
suppressPackageStartupMessages({ library(ragg); library(png) })
FONT <- { pref <- c("Helvetica","Arial","Liberation Sans","DejaVu Sans")
  fams <- tryCatch(systemfonts::system_fonts()$family, error=function(e) character(0))
  hit <- pref[pref %in% fams]; if (length(hit)) hit[1] else "sans" }
.trimS <- function(img, tol = 0.985, pad = 5) {
  d <- dim(img); if (length(d) == 2) img <- array(img, c(d, 1))
  ch <- min(dim(img)[3], 3); ink <- Reduce(`|`, lapply(seq_len(ch), function(k) img[, , k] < tol))
  r <- which(rowSums(ink) > 0); c <- which(colSums(ink) > 0)
  if (!length(r) || !length(c)) return(img)
  img[max(1,min(r)-pad):min(dim(img)[1],max(r)+pad),
      max(1,min(c)-pad):min(dim(img)[2],max(c)+pad), , drop = FALSE]
}
p_bare <- fig_activity_bar(D$act, n = 14, mode = "publication") +
  theme_publication() + theme(text = element_text(family = FONT), plot.title = element_blank(),
    axis.title = element_text(size = 14), axis.text = element_text(size = 11.5),
    legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
    plot.margin = margin(3, 4, 3, 4))
dir.create("paper/panels/figureS1", showWarnings = FALSE, recursive = TRUE)
f <- "paper/panels/figureS1/a_activity.png"
ggsave(f, p_bare, width = 7.5, height = 5.0, dpi = 400, bg = "white", device = ragg::agg_png)
png::writePNG(.trimS(png::readPNG(f)), f)
cat("wrote figureS1 bare panel\n")

## ---- report the top pathways so the caption/text can quote real numbers
top <- head(D$act[order(-abs(D$act$score)), ], 6)
cat("\nTop pathways by |activity score| (fill these into paper.md TODO):\n")
print(top[, intersect(c("source", "score", "p_value", "padj"), colnames(top))],
      row.names = FALSE)
