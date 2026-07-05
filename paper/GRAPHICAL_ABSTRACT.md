# Figure 1 — Graphical abstract (BioRender-style)

The manuscript's **Figure 1** slot is an overview / graphical abstract. Build it
in BioRender (or an image generator) and drop the export in as
`paper/figures/figure1.png` — the caption ("Overview of RNAflow") is already
written. Target ~180 mm wide, landscape, flat scientific style.

## The one-line message
**RNAflow: one tool, from a count matrix to a reproducible, multi-layered
interpretation.**

## Must contain (left → right flow)

1. **Input (left).** A count-matrix / spreadsheet icon labelled "Counts +
   metadata", with three small organism icons (human, mouse, rat).
2. **Analyses (centre).** A hub / pipeline of the **nine modules**, each with a
   tiny representative glyph:
   - Differential expression — volcano plot
   - Quality control — MA plot / gauge
   - Sample overview — PCA scatter (PCA / UMAP / 3D)
   - Functional enrichment — dot/bar (GSEA + ORA)
   - Co-expression network — WGCNA graph
   - Regulator & pathway activity — TF → gene arrows (decoupleR)
   - Per-sample signatures — heatmap (GSVA)
   - Multi-contrast comparison — Venn / UpSet
   - AI-assisted interpretation — chat/AI chip (optional)
   Group visually into 3 bands (exploration · interpretation · synthesis).
3. **Outputs (right).** "Reproducible export": R-script icon, a document
   (Methods paragraph), a report/page, a session file (.rds), and a Docker whale.

Arrows: Input → Analyses → Outputs. A subtle loop arrow from Outputs back to
Input conveys "reproducible / re-runnable".

## Style
- Flat, clean, BioRender scientific look; rounded cards; soft palette.
- Accent green **#1D9E75** (RNAflow), secondary blue **#4A90D9**, warm
  **#D55E00**; neutral greys for boxes; white background.
- Sans-serif (Inter / Helvetica). No 3D, no gradients-heavy, no clutter.
- Colour-vision-safe; readable at 1-column width.

## Ready-to-use prompt (adapt as you like)

> A clean, flat scientific **graphical abstract** for a bioinformatics software
> tool called **RNAflow**, landscape, white background, BioRender style, soft
> professional palette (green #1D9E75, blue #4A90D9, orange #D55E00, greys).
> Left: an icon of a gene count matrix / spreadsheet labelled "Counts +
> metadata" with small human, mouse and rat icons. A big arrow to the centre.
> Centre: a titled panel "Analyses" containing nine small rounded cards, each
> with a tiny icon and label — "Differential expression" (volcano plot), "Quality
> control" (MA plot), "Sample overview PCA/UMAP" (scatter clusters), "Functional
> enrichment" (dot plot), "Co-expression WGCNA" (network graph), "Activity
> decoupleR" (transcription-factor to gene arrows), "Per-sample signatures GSVA"
> (heatmap), "Multi-contrast comparison" (Venn/UpSet), "AI interpretation" (AI
> chip). A big arrow to the right. Right: a panel "Reproducible export" with
> icons for an R script, a Methods paragraph document, an HTML report, a session
> file, and a Docker whale. A thin curved arrow loops from the outputs back to
> the input to suggest reproducibility. Minimal, uncluttered, publication-ready,
> readable labels, no photorealism.

## Slot it in
Save the final artwork as **`paper/figures/figure1.png`** (and ideally a vector
`figure1.pdf`). The placeholder schematic currently at that path (from
`paper/make_overview.R`) is only a stand-in.
