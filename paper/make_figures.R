# Generate Figure 1 for the RNAflow manuscript from the two bundled demo
# datasets, and print the exact statistics quoted in the Results.
# Run from the package root:  Rscript paper/make_figures.R
suppressPackageStartupMessages({
  if (!suppressWarnings(require(RNAflow, quietly = TRUE))) {
    devtools::load_all(".", quiet = TRUE)          # dev: package not installed
  }
  library(ggplot2); library(patchwork)
})
set.seed(1)

# Okabe-Ito 8-colour palette (colour-vision-deficiency safe, fixed order)
OI <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
        "#0072B2", "#D55E00", "#CC79A7", "#7F7F7F")

ext <- function(f) system.file("extdata", f, package = "RNAflow")

## ---- airway (simple 2-group) -------------------------------------------
a_counts <- read_counts(ext("demo_airway_counts.csv"))
a_meta   <- read_metadata(ext("demo_airway_metadata.csv"),
                          counts_samples = colnames(a_counts))
a_res <- run_deseq2(a_counts, a_meta, design = ~ cell + condition,
                    contrast = c("condition", "Dex", "Control"), shrink = TRUE)
a_sig <- a_res[!is.na(a_res$padj) & a_res$padj < 0.05 & abs(a_res$log2FoldChange) > 1, ]
cat(sprintf("[airway] %d genes tested; %d significant (padj<0.05,|LFC|>1); up=%d down=%d\n",
            nrow(a_res), nrow(a_sig),
            sum(a_sig$log2FoldChange > 0), sum(a_sig$log2FoldChange < 0)))
top_up <- utils::head(a_sig$gene[order(-a_sig$log2FoldChange)], 6)
cat("[airway] top up:", paste(top_up, collapse = ", "), "\n")

pA <- fig_volcano(a_res, lfc_thr = 1, padj_thr = 0.05, n_label = 12,
                  mode = "publication") +
  ggtitle("airway: dexamethasone vs control")

## ---- TCGA (complex 8-group) --------------------------------------------
t_counts <- read_counts(ext("demo_tcga_counts.csv"))
t_meta   <- read_metadata(ext("demo_tcga_metadata.csv"),
                          counts_samples = colnames(t_counts))
vst <- normalize_counts(t_counts, t_meta, "vst")

pca <- compute_pca(vst, 500)
sc  <- merge(pca$scores, t_meta, by = "sample")
cat(sprintf("[TCGA] PCA PC1=%.1f%% PC2=%.1f%%\n", pca$pct[1], pca$pct[2]))
pB <- ggplot(sc, aes(PC1, PC2, colour = cancer_type)) +
  geom_point(size = 2.4, alpha = 0.9) +
  scale_colour_manual(values = OI, name = "Cancer type") +
  labs(x = sprintf("PC1 (%.1f%%)", pca$pct[1]),
       y = sprintf("PC2 (%.1f%%)", pca$pct[2]),
       title = "TCGA pan-cancer: PCA") +
  theme_publication()

um <- compute_umap(vst, 500, seed = 42)
uc <- merge(um$scores, t_meta, by = "sample")
xy <- grep("UMAP", colnames(uc), value = TRUE)
pC <- ggplot(uc, aes(.data[[xy[1]]], .data[[xy[2]]], colour = cancer_type)) +
  geom_point(size = 2.4, alpha = 0.9) +
  scale_colour_manual(values = OI, name = "Cancer type") +
  labs(x = "UMAP1", y = "UMAP2", title = "TCGA pan-cancer: UMAP") +
  theme_publication()

# one representative pairwise contrast (glioma vs lung adenocarcinoma)
t_res <- run_deseq2(t_counts, t_meta, design = ~ cancer_type,
                    contrast = c("cancer_type", "LGG", "LUAD"), shrink = TRUE)
t_sig <- t_res[!is.na(t_res$padj) & t_res$padj < 0.05 & abs(t_res$log2FoldChange) > 1, ]
cat(sprintf("[TCGA] LGG vs LUAD: %d significant (padj<0.05,|LFC|>1)\n", nrow(t_sig)))

## ---- functional enrichment panel (airway GSEA, Hallmark) ----------------
pD <- tryCatch({
  sets <- get_gene_sets("human", collection = "H")
  gsea <- run_gsea(a_res, sets, rank_by = "stat")
  cat(sprintf("[airway] GSEA Hallmark: %d sets at padj<0.05\n",
              sum(gsea$padj < 0.05, na.rm = TRUE)))
  fig_enrich_dot(gsea, n = 10) + ggtitle("airway: GSEA (MSigDB Hallmark)")
}, error = function(e) {
  cat("[enrichment panel skipped]:", conditionMessage(e), "\n")
  fig_volcano(t_res, lfc_thr = 1, padj_thr = 0.05, n_label = 10,
              mode = "publication") + ggtitle("TCGA: LGG vs LUAD")
})

## ---- compose -----------------------------------------------------------
fig <- (pB | pC) / (pA | pD) +
  patchwork::plot_annotation(tag_levels = "A")
ggsave("paper/figures/figure1.png", fig, width = 12, height = 9.5, dpi = 300,
       bg = "white")
ggsave("paper/figures/figure1.pdf", fig, width = 12, height = 9.5, bg = "white")
cat("Wrote paper/figures/figure1.{png,pdf}\n")
