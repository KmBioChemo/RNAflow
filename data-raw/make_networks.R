# data-raw/make_networks.R
#
# Regenerate the offline prior-knowledge networks bundled with RNAflow:
#
#   inst/extdata/collectri_human.rds   CollecTRI TF regulons  (source, target, mor)
#   inst/extdata/progeny_human.rds     PROGENy pathway footprints
#                                      (source, target, weight, p_value)
#
# Why bundle them: activity inference (R/analysis_decoupler.R) normally asks
# decoupleR to download these from the OmniPath web service via OmnipathR.
# That path is fragile -- the service is sometimes down and its offline
# static-table fallback is broken upstream -- so CollecTRI activity often
# never worked for users. get_tf_network() / get_pathway_network() now fall
# back to these bundled copies (human only) when the live fetch fails.
#
# These snapshots are downloaded directly from the OmniPath REST API rather
# than through decoupleR/OmnipathR, so this script has NO Bioconductor
# dependency -- only base R (utils::download.file, jsonlite is not required
# because we request TSV). Run it from the package root:
#
#   Rscript data-raw/make_networks.R
#
# Only human (organism 9606) is snapshotted; other organisms still require the
# live OmniPath fetch. Re-run to refresh against the current OmniPath release.

dir.create("inst/extdata", showWarnings = FALSE, recursive = TRUE)

read_omnipath_tsv <- function(url) {
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  utils::download.file(url, tmp, quiet = TRUE)
  utils::read.delim(tmp, sep = "\t", header = TRUE,
                    stringsAsFactors = FALSE, check.names = FALSE)
}

## ---- CollecTRI (transcription-factor regulons) ---------------------------
# decoupleR's mode-of-regulation (mor) is +1 for activation and -1 for pure
# repression. OmniPath exposes is_stimulation / is_inhibition per interaction.
message("Downloading CollecTRI ...")
collectri_url <- paste0(
  "https://omnipathdb.org/interactions",
  "?datasets=collectri",
  "&fields=sources,references",
  "&genesymbols=1",
  "&organisms=9606",
  "&license=academic")

ct <- read_omnipath_tsv(collectri_url)

# Map to decoupleR's (source, target, mor) schema.
ct$mor <- ifelse(ct$is_inhibition == "True" & ct$is_stimulation != "True",
                 -1L, 1L)
collectri <- data.frame(
  source = ct$source_genesymbol,
  target = ct$target_genesymbol,
  mor    = ct$mor,
  stringsAsFactors = FALSE)

# Drop blanks and duplicate edges (keep first).
collectri <- collectri[nzchar(collectri$source) & nzchar(collectri$target), ]
collectri <- collectri[!duplicated(collectri[c("source", "target")]), ]
rownames(collectri) <- NULL

saveRDS(collectri, "inst/extdata/collectri_human.rds")
message(sprintf("  wrote collectri_human.rds: %d edges, %d TFs, %d repressions",
                nrow(collectri), length(unique(collectri$source)),
                sum(collectri$mor < 0)))

## ---- PROGENy (pathway-responsive genes) ----------------------------------
# PROGENy is served as long-format annotations (one row per record_id/label);
# pivot back to one row per (pathway, gene) with its weight and p_value, then
# keep the top-500 most responsive genes per pathway (smallest p_value), which
# is decoupleR::get_progeny(top = 500)'s default footprint size.
message("Downloading PROGENy ...")
progeny_url <- paste0(
  "https://omnipathdb.org/annotations",
  "?resources=PROGENy",
  "&license=academic")

pr <- read_omnipath_tsv(progeny_url)

# Long -> wide by record_id: each record has pathway / weight / p_value / genesymbol.
wide <- reshape(
  pr[, c("record_id", "label", "value")],
  idvar = "record_id", timevar = "label", direction = "wide")
colnames(wide) <- sub("^value\\.", "", colnames(wide))

progeny <- data.frame(
  source  = wide$pathway,
  target  = wide$genesymbol,
  weight  = as.numeric(wide$weight),
  p_value = as.numeric(wide$p_value),
  stringsAsFactors = FALSE)
progeny <- progeny[nzchar(progeny$source) & nzchar(progeny$target) &
                   !is.na(progeny$weight) & !is.na(progeny$p_value), ]

# Top 500 per pathway by ascending p_value.
progeny <- do.call(rbind, lapply(split(progeny, progeny$source), function(d) {
  d <- d[order(d$p_value), , drop = FALSE]
  utils::head(d, 500)
}))
rownames(progeny) <- NULL

saveRDS(progeny, "inst/extdata/progeny_human.rds")
message(sprintf("  wrote progeny_human.rds: %d edges, %d pathways",
                nrow(progeny), length(unique(progeny$source))))

message("Done.")
