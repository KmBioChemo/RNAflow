# Validation analyses for the RNAflow manuscript, and Figure 4:
#   1. Reproducibility round-trip: export the R script, re-run it in a clean R
#      process, and check the DE table is identical.
#   2. All-pairwise consistency: contrasts from one shared fit == independent fits.
#   3. Method concordance: RNAflow (DESeq2) vs limma-voom on airway.
# (Test coverage is measured separately with covr.)
# Run from the package root:  Rscript paper/make_validation.R
suppressPackageStartupMessages({
  if (!suppressWarnings(require(RNAflow, quietly = TRUE))) devtools::load_all(".", quiet = TRUE)
  library(ggplot2); library(patchwork); library(ragg)
})
set.seed(1)
FONT <- { pref <- c("Helvetica","Arial","Liberation Sans","DejaVu Sans")
  fams <- tryCatch(systemfonts::system_fonts()$family, error=function(e) character(0))
  hit <- pref[pref %in% fams]; if (length(hit)) hit[1] else "sans" }
ext <- function(f) system.file("extdata", f, package = "RNAflow")
idline <- geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "#B0B7C0")

ac <- read_counts(ext("demo_airway_counts.csv"))
am <- read_metadata(ext("demo_airway_metadata.csv"), counts_samples = colnames(ac))
res_orig <- run_deseq2(ac, am, design = ~ cell + condition,
                       contrast = c("condition", "Dex", "Control"), shrink = TRUE)

## ===== 1. Reproducibility round-trip ===================================
p <- empty_project("airway")
p$counts <- ac; p$metadata <- am; p$organism <- "human"
p$contrasts <- contrast_store_upsert(
  list(), "condition: Dex vs Control", res_orig,
  params = list(design_var = "condition", treated = "Dex", reference = "Control",
                covariates = "cell", shrink = TRUE))
script <- generate_r_script(p, counts_path = ext("demo_airway_counts.csv"),
                            metadata_path = ext("demo_airway_metadata.csv"))
script <- unlist(strsplit(paste(script, collapse = "\n"), "\n"))
script <- sub("^library\\(RNAflow\\)$",
              "suppressPackageStartupMessages(devtools::load_all('.', quiet=TRUE))", script)
out_csv <- file.path(tempdir(), "roundtrip_out.csv")
tmpR <- tempfile(fileext = ".R")
writeLines(c(script,
             'robj <- mget(ls(pattern="^res_"))[[1]]',
             sprintf('write.csv(robj, "%s", row.names=FALSE)', out_csv)), tmpR)
rc <- system2(file.path(R.home("bin"), "Rscript"), tmpR, stdout = FALSE, stderr = FALSE)
res_rep <- read.csv(out_csv, stringsAsFactors = FALSE)
m1 <- merge(res_orig[, c("gene", "log2FoldChange", "padj")],
            res_rep[, c("gene", "log2FoldChange", "padj")], by = "gene",
            suffixes = c(".o", ".r"))
rt_maxdiff <- max(abs(m1$log2FoldChange.o - m1$log2FoldChange.r), na.rm = TRUE)
rt_cor <- stats::cor(m1$log2FoldChange.o, m1$log2FoldChange.r, use = "complete.obs")
cat(sprintf("[ROUNDTRIP] genes=%d  max|dLFC|=%.2e  r=%.6f  identical=%s\n",
            nrow(m1), rt_maxdiff, rt_cor, rt_maxdiff < 1e-8))

## ===== 2. All-pairwise consistency =====================================
tc <- read_counts(ext("demo_tcga_counts.csv"))
tm <- read_metadata(ext("demo_tcga_metadata.csv"), counts_samples = colnames(tc))
lst <- run_deseq2_all_pairs(tc, tm, design = ~ cancer_type, shrink = FALSE)
sel <- c("cancer_type: LUAD vs LGG", "cancer_type: COAD vs BRCA", "cancer_type: THCA vs KIRC")
sel <- intersect(sel, names(lst)); if (length(sel) < 3) sel <- names(lst)[1:3]
maxd <- 0; ap_df <- NULL
for (lab in sel) {
  sh <- lst[[lab]]; pr <- attr(sh, "pair")     # c(treated, reference)
  ind <- run_deseq2(tc, tm, design = ~ cancer_type,
                    contrast = c("cancer_type", pr["treated"], pr["reference"]), shrink = FALSE)
  mm <- merge(sh[, c("gene","log2FoldChange","stat")],
              ind[, c("gene","log2FoldChange","stat")], by = "gene", suffixes = c(".s",".i"))
  maxd <- max(maxd, max(abs(mm$log2FoldChange.s - mm$log2FoldChange.i), na.rm = TRUE),
              max(abs(mm$stat.s - mm$stat.i), na.rm = TRUE))
  if (is.null(ap_df)) ap_df <- mm
}
ap_cor <- stats::cor(ap_df$log2FoldChange.s, ap_df$log2FoldChange.i, use = "complete.obs")
cat(sprintf("[ALLPAIRS] %d pairs  max|d(LFC,stat)|=%.2e  r(first)=%.6f\n",
            length(sel), maxd, ap_cor))

## ===== 3. Concordance with limma-voom (airway) =========================
suppressPackageStartupMessages(library(limma))
keep <- rowSums(ac >= 10) >= 4
cd <- data.frame(cell = factor(am$cell), condition = factor(am$condition, c("Control","Dex")))
design <- stats::model.matrix(~ cell + condition, data = cd)
v <- limma::voom(ac[keep, ], design)
fit <- limma::eBayes(limma::lmFit(v, design))
tt <- limma::topTable(fit, coef = "conditionDex", number = Inf, sort.by = "none")
tt$gene <- rownames(tt)
mc <- merge(res_orig[, c("gene","log2FoldChange","padj")],
            tt[, c("gene","logFC","adj.P.Val")], by = "gene")
co_r  <- stats::cor(mc$log2FoldChange, mc$logFC, use = "complete.obs")
co_rs <- stats::cor(mc$log2FoldChange, mc$logFC, method = "spearman", use = "complete.obs")
sd <- mc$gene[!is.na(mc$padj) & mc$padj < 0.05 & abs(mc$log2FoldChange) > 1]
sl <- mc$gene[!is.na(mc$adj.P.Val) & mc$adj.P.Val < 0.05 & abs(mc$logFC) > 1]
jac <- length(intersect(sd, sl)) / length(union(sd, sl))
cat(sprintf("[CONCORDANCE] n=%d  Pearson=%.3f  Spearman=%.3f  Jaccard(sig)=%.3f (DE:%d limma:%d shared:%d)\n",
            nrow(mc), co_r, co_rs, jac, length(sd), length(sl), length(intersect(sd, sl))))

## ===== Figure 5 bare panels for the figure pipeline =========================
# The composed Figure 5 is assembled by paper/figure5_rebuild/ from these bare
# panels; this script only computes the validation and exports the panels.
comma <- function(n) formatC(n, big.mark = ",", format = "d")
# Rebuild the three scatters with no titles (Python owns lettering/titles),
# larger axis fonts, and export tight-cropped like the other figures' panels.
theme_bare5 <- theme_publication() + theme(
  text = element_text(family = FONT),
  plot.title = element_blank(), plot.subtitle = element_blank(),
  axis.title = element_text(size = 14), axis.text = element_text(size = 11.5),
  plot.margin = margin(3, 4, 3, 4))
annot_b <- function(lbl) annotate("label", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.06,
  label = lbl, size = 4.4, label.size = 0, fill = "#ffffffcc", colour = "#2b2f36", lineheight = 0.98)
.trim5 <- function(img, tol = 0.985, pad = 5) {
  d <- dim(img); if (length(d) == 2) img <- array(img, c(d, 1))
  ch <- min(dim(img)[3], 3); ink <- Reduce(`|`, lapply(seq_len(ch), function(k) img[, , k] < tol))
  r <- which(rowSums(ink) > 0); c <- which(colSums(ink) > 0)
  if (!length(r) || !length(c)) return(img)
  img[max(1,min(r)-pad):min(dim(img)[1],max(r)+pad),
      max(1,min(c)-pad):min(dim(img)[2],max(c)+pad), , drop = FALSE]
}
save_bare5 <- function(p, name, w = 4.9, h = 4.8) {
  dir.create("paper/panels/figure5", showWarnings = FALSE, recursive = TRUE)
  f <- file.path("paper/panels/figure5", paste0(name, ".png"))
  ggsave(f, p, width = w, height = h, dpi = 400, bg = "white", device = ragg::agg_png)
  png::writePNG(.trim5(png::readPNG(f)), f); cat("  fig5", name, "\n")
}
pA_b <- ggplot(m1, aes(log2FoldChange.o, log2FoldChange.r)) + idline +
  geom_point(size = 1.2, alpha = .4, colour = "#1D9E75") +
  annot_b(sprintf("r = %.4f\nmax abs diff = %.0e\nn = %s genes", rt_cor, rt_maxdiff, comma(nrow(m1)))) +
  labs(x = "log2FC (interactive)", y = "log2FC (re-run script)") + theme_bare5
pB_b <- ggplot(ap_df, aes(log2FoldChange.i, log2FoldChange.s)) + idline +
  geom_point(size = 1.2, alpha = .4, colour = "#0072B2") +
  annot_b(sprintf("r = %.4f\nmax abs diff = %.0f\n3 contrasts", ap_cor, maxd)) +
  labs(x = "log2FC (independent fit)", y = "log2FC (shared fit)") + theme_bare5
pC_b <- ggplot(mc, aes(logFC, log2FoldChange)) + idline +
  geom_point(size = 1.2, alpha = .4, colour = "#D55E00") +
  annot_b(sprintf("Pearson r = %.3f\nSpearman = %.3f\nJaccard(sig) = %.2f\nn = %s genes",
                  co_r, co_rs, jac, comma(nrow(mc)))) +
  labs(x = "log2FC (limma-voom)", y = "log2FC (RNAflow / DESeq2)") + theme_bare5
suppressPackageStartupMessages(library(png))
save_bare5(pA_b, "a_roundtrip"); save_bare5(pB_b, "b_allpairwise"); save_bare5(pC_b, "c_concordance")
cat("Wrote figure5 bare panels\n")
