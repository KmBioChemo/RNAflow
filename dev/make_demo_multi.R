# Generate a richer, factorial demo dataset for RNAflow.
#
# Design: genotype {WT, KO} x treatment {Veh, LPS, LPS_Drug}, n = 4 per cell
#   -> 24 samples, plus a nuisance `batch` covariate.
# A `group` column (6 levels) lets the single-term DE module run ANY pairwise
# contrast, so the multi-contrast Compare tab has plenty to chew on.
#
# Planted signal (controlled overlaps for Venn / UpSet):
#   M_infl   : up with LPS in both genotypes (inflammation)
#   M_rescue : up with LPS, brought back down by the drug (overlaps "LPS vs Veh"
#              and "LPS_Drug vs LPS")
#   M_geno   : KO vs WT shift, treatment-independent
#   M_inter  : LPS effect present ONLY in KO (genotype x treatment interaction)
#
# Reproducible: fixed seed. Writes CSVs to inst/extdata/.

set.seed(2026)

n_genes <- 4000

## ---- gene universe: real mouse symbols when available -----------------------
symbols <- NULL
if (requireNamespace("org.Mm.eg.db", quietly = TRUE) &&
    requireNamespace("AnnotationDbi", quietly = TRUE)) {
  all_sym <- AnnotationDbi::keys(org.Mm.eg.db::org.Mm.eg.db, keytype = "SYMBOL")
  all_sym <- unique(all_sym[grepl("^[A-Za-z][A-Za-z0-9]+$", all_sym)])
  if (length(all_sym) >= n_genes) {
    symbols <- sample(all_sym, n_genes)
    message("Using ", n_genes, " real mouse symbols from org.Mm.eg.db")
  }
}
if (is.null(symbols)) {
  symbols <- sprintf("Gene%04d", seq_len(n_genes))
  message("org.Mm.eg.db not available - using synthetic gene IDs")
}

## ---- experimental design ----------------------------------------------------
genotypes  <- c("WT", "KO")
treatments <- c("Veh", "LPS", "LPS_Drug")
n_rep      <- 4

design <- expand.grid(rep = seq_len(n_rep),
                      treatment = treatments,
                      genotype  = genotypes,
                      stringsAsFactors = FALSE)
design$group  <- paste(design$genotype, design$treatment, sep = "_")
design$sample <- sprintf("%s_%s_%d", design$genotype, design$treatment, design$rep)
# Balanced nuisance batch (confounded with neither factor)
design$batch  <- rep(c("B1", "B2"), length.out = nrow(design))
n_samp <- nrow(design)

groups <- c("WT_Veh", "WT_LPS", "WT_LPS_Drug", "KO_Veh", "KO_LPS", "KO_LPS_Drug")

## ---- planted log2 fold changes (genes x groups), baseline = WT_Veh ----------
lfc <- matrix(0, nrow = n_genes, ncol = length(groups),
              dimnames = list(symbols, groups))

idx <- sample(seq_len(n_genes))
take <- function(n) { picked <- idx[seq_len(n)]; idx <<- idx[-seq_len(n)]; picked }

m_infl   <- take(300)   # inflammation, both genotypes
m_rescue <- take(150)   # LPS up, drug rescues
m_geno   <- take(200)   # genotype effect
m_inter  <- take(120)   # LPS effect only in KO

set_lfc <- function(rows, vals) for (g in names(vals)) lfc[rows, g] <<- vals[[g]]

# Inflammation: up with LPS (drug keeps most of it)
set_lfc(m_infl, list(
  WT_LPS = rnorm(length(m_infl),  2.0, 0.4),
  WT_LPS_Drug = rnorm(length(m_infl), 1.7, 0.4),
  KO_LPS = rnorm(length(m_infl),  2.1, 0.4),
  KO_LPS_Drug = rnorm(length(m_infl), 1.8, 0.4)))

# Rescue: up with LPS, returned to baseline by the drug
set_lfc(m_rescue, list(
  WT_LPS = rnorm(length(m_rescue),  2.2, 0.4),
  WT_LPS_Drug = rnorm(length(m_rescue), 0.1, 0.3),
  KO_LPS = rnorm(length(m_rescue),  2.0, 0.4),
  KO_LPS_Drug = rnorm(length(m_rescue), 0.8, 0.4)))

# Genotype: KO shifted vs WT across all treatments
geno_shift <- rnorm(length(m_geno), 1.6, 0.5) * sample(c(1, -1), length(m_geno), TRUE)
set_lfc(m_geno, list(
  KO_Veh = geno_shift, KO_LPS = geno_shift, KO_LPS_Drug = geno_shift))

# Interaction: LPS induces these only in KO
set_lfc(m_inter, list(
  KO_LPS = rnorm(length(m_inter), 2.3, 0.4),
  KO_LPS_Drug = rnorm(length(m_inter), 2.1, 0.4)))

## ---- simulate counts (negative binomial) ------------------------------------
# Gene baseline expression (broad, heavy-tailed) and per-gene dispersion.
base_mu   <- 2^(rnorm(n_genes, 5.5, 2.2))            # baseline mean in WT_Veh
disp      <- 0.10 + 3 / (base_mu + 8)               # higher dispersion for low counts
size      <- 1 / disp                               # NB size param
lib_factor <- exp(rnorm(n_samp, 0, 0.18))           # library-size variation
gene_batch_shift <- rnorm(n_genes, 0, 0.25)         # mild per-gene B2 shift (log2)

counts <- matrix(0L, nrow = n_genes, ncol = n_samp,
                 dimnames = list(symbols, design$sample))
for (j in seq_len(n_samp)) {
  g <- design$group[j]
  batch_mult <- if (design$batch[j] == "B2") 2^gene_batch_shift else rep(1, n_genes)
  mu <- base_mu * 2^(lfc[, g]) * lib_factor[j] * batch_mult
  counts[, j] <- rnbinom(n_genes, mu = pmax(mu, 1e-3), size = size)
}

## ---- write ------------------------------------------------------------------
counts_df <- data.frame(gene_id = rownames(counts), counts, check.names = FALSE)
meta_df   <- design[, c("sample", "group", "genotype", "treatment", "batch")]

out_dir <- file.path("inst", "extdata")
write.csv(counts_df, file.path(out_dir, "demo_multi_counts.csv"), row.names = FALSE)
write.csv(meta_df,   file.path(out_dir, "demo_multi_metadata.csv"), row.names = FALSE)

message(sprintf("Wrote %d genes x %d samples (%d groups, batch covariate).",
                n_genes, n_samp, length(groups)))
message("Planted modules: infl=", length(m_infl), " rescue=", length(m_rescue),
        " geno=", length(m_geno), " inter=", length(m_inter))
