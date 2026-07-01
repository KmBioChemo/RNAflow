# Package index

## Application

Launch the interactive Shiny app.

- [`run_app()`](https://KmBioChemo.github.io/RNAflow/reference/run_app.md)
  : Launch the RNAflow Shiny application

- [`launch_app()`](https://KmBioChemo.github.io/RNAflow/reference/launch_app.md)
  :

  Alias for
  [`run_app()`](https://KmBioChemo.github.io/RNAflow/reference/run_app.md)

## Data input & validation

Read and strictly validate counts, metadata, and DE results.

- [`read_counts()`](https://KmBioChemo.github.io/RNAflow/reference/read_counts.md)
  : Read a counts matrix from a file
- [`read_metadata()`](https://KmBioChemo.github.io/RNAflow/reference/read_metadata.md)
  : Read sample metadata from a file
- [`read_de_results()`](https://KmBioChemo.github.io/RNAflow/reference/read_de_results.md)
  : Read a DE results table from a file
- [`validate_counts()`](https://KmBioChemo.github.io/RNAflow/reference/validate_counts.md)
  : Validate a counts matrix
- [`validate_metadata()`](https://KmBioChemo.github.io/RNAflow/reference/validate_metadata.md)
  : Validate sample metadata
- [`validate_de_results()`](https://KmBioChemo.github.io/RNAflow/reference/validate_de_results.md)
  : Validate a DE results table

## Differential expression

- [`run_deseq2()`](https://KmBioChemo.github.io/RNAflow/reference/run_deseq2.md)
  : Run DESeq2 on a counts matrix
- [`normalize_counts()`](https://KmBioChemo.github.io/RNAflow/reference/normalize_counts.md)
  : Normalize counts (variance-stabilized transform)

## Core figures

Volcano, heatmap and PCA, with exploration / publication modes.

- [`fig_volcano()`](https://KmBioChemo.github.io/RNAflow/reference/fig_volcano.md)
  : Volcano plot functions
- [`fig_volcano_interactive()`](https://KmBioChemo.github.io/RNAflow/reference/fig_volcano_interactive.md)
  : Interactive volcano plot (plotly)
- [`fig_heatmap()`](https://KmBioChemo.github.io/RNAflow/reference/fig_heatmap.md)
  : Heatmap figures
- [`fig_pca()`](https://KmBioChemo.github.io/RNAflow/reference/fig_pca.md)
  : PCA figures
- [`compute_pca()`](https://KmBioChemo.github.io/RNAflow/reference/compute_pca.md)
  : Compute PCA scores

## QC & diagnostics

- [`fig_pval_hist()`](https://KmBioChemo.github.io/RNAflow/reference/fig_pval_hist.md)
  : P-value histogram
- [`fig_ma()`](https://KmBioChemo.github.io/RNAflow/reference/fig_ma.md)
  : MA plot
- [`fig_sample_cor()`](https://KmBioChemo.github.io/RNAflow/reference/fig_sample_cor.md)
  : Sample-sample correlation heatmap
- [`fig_lib_sizes()`](https://KmBioChemo.github.io/RNAflow/reference/fig_lib_sizes.md)
  : Library-size bar chart

## Themes & figure export

- [`fig_theme()`](https://KmBioChemo.github.io/RNAflow/reference/fig_theme.md)
  : Figure theme system
- [`theme_publication()`](https://KmBioChemo.github.io/RNAflow/reference/theme_publication.md)
  : Publication-ready ggplot2 theme
- [`theme_exploration()`](https://KmBioChemo.github.io/RNAflow/reference/theme_exploration.md)
  : Exploration-mode ggplot2 theme
- [`save_ggplot()`](https://KmBioChemo.github.io/RNAflow/reference/save_ggplot.md)
  : Save a ggplot to disk
- [`save_pheatmap()`](https://KmBioChemo.github.io/RNAflow/reference/save_pheatmap.md)
  : Save a pheatmap to disk
- [`save_compare()`](https://KmBioChemo.github.io/RNAflow/reference/save_compare.md)
  : Save a comparison figure to disk

## Multi-contrast comparison

Compare significant-gene sets and signatures across contrasts.

- [`contrast_sig_genes()`](https://KmBioChemo.github.io/RNAflow/reference/contrast_sig_genes.md)
  : Significant genes of a single contrast
- [`contrast_sig_sets()`](https://KmBioChemo.github.io/RNAflow/reference/contrast_sig_sets.md)
  : Significant-gene sets across contrasts
- [`contrast_lfc_matrix()`](https://KmBioChemo.github.io/RNAflow/reference/contrast_lfc_matrix.md)
  : Gene x contrast log2FoldChange matrix
- [`fig_venn()`](https://KmBioChemo.github.io/RNAflow/reference/fig_venn.md)
  : Venn diagram of significant-gene sets
- [`fig_upset()`](https://KmBioChemo.github.io/RNAflow/reference/fig_upset.md)
  : UpSet plot of significant-gene sets
- [`fig_volcano_grid()`](https://KmBioChemo.github.io/RNAflow/reference/fig_volcano_grid.md)
  : Side-by-side volcano grid
- [`fig_lfc_heatmap()`](https://KmBioChemo.github.io/RNAflow/reference/fig_lfc_heatmap.md)
  : Cross-contrast log2FoldChange signature heatmap

## Functional enrichment

GSEA (fgsea / MSigDB) and ORA (clusterProfiler; GO / KEGG / Reactome).

- [`get_gene_sets()`](https://KmBioChemo.github.io/RNAflow/reference/get_gene_sets.md)
  : Fetch MSigDB gene sets for an organism
- [`rank_genes()`](https://KmBioChemo.github.io/RNAflow/reference/rank_genes.md)
  : Build a ranked gene vector from DE results
- [`run_gsea()`](https://KmBioChemo.github.io/RNAflow/reference/run_gsea.md)
  : Run GSEA (fgsea) on DE results
- [`run_ora()`](https://KmBioChemo.github.io/RNAflow/reference/run_ora.md)
  : Run over-representation analysis (ORA)
- [`fig_enrich_dot()`](https://KmBioChemo.github.io/RNAflow/reference/fig_enrich_dot.md)
  : Enrichment dotplot
- [`fig_enrich_bar()`](https://KmBioChemo.github.io/RNAflow/reference/fig_enrich_bar.md)
  : Enrichment lollipop bar
- [`fig_gsea_curve()`](https://KmBioChemo.github.io/RNAflow/reference/fig_gsea_curve.md)
  : GSEA running-enrichment curve
- [`fig_enrich_map()`](https://KmBioChemo.github.io/RNAflow/reference/fig_enrich_map.md)
  : Enrichment map (pathway network)
- [`fig_gsea_ridge()`](https://KmBioChemo.github.io/RNAflow/reference/fig_gsea_ridge.md)
  : GSEA ridgeline plot

## Co-expression networks (WGCNA)

- [`wgcna_datexpr()`](https://KmBioChemo.github.io/RNAflow/reference/wgcna_datexpr.md)
  : Build the WGCNA expression matrix
- [`wgcna_pick_power()`](https://KmBioChemo.github.io/RNAflow/reference/wgcna_pick_power.md)
  : Soft-threshold (power) selection
- [`run_wgcna()`](https://KmBioChemo.github.io/RNAflow/reference/run_wgcna.md)
  : Detect co-expression modules
- [`build_traits()`](https://KmBioChemo.github.io/RNAflow/reference/build_traits.md)
  : Build a numeric trait matrix from sample metadata
- [`module_trait_cor()`](https://KmBioChemo.github.io/RNAflow/reference/module_trait_cor.md)
  : Module-trait correlation
- [`hub_genes()`](https://KmBioChemo.github.io/RNAflow/reference/hub_genes.md)
  : Intramodular hub genes
- [`module_gene_list()`](https://KmBioChemo.github.io/RNAflow/reference/module_gene_list.md)
  : Module gene lists
- [`module_summary()`](https://KmBioChemo.github.io/RNAflow/reference/module_summary.md)
  : Per-module summary table
- [`fig_soft_threshold()`](https://KmBioChemo.github.io/RNAflow/reference/fig_soft_threshold.md)
  : Soft-threshold diagnostic plot
- [`fig_module_trait()`](https://KmBioChemo.github.io/RNAflow/reference/fig_module_trait.md)
  : Module-trait correlation heatmap
- [`fig_module_sizes()`](https://KmBioChemo.github.io/RNAflow/reference/fig_module_sizes.md)
  : Module size bar chart
- [`fig_eigengene()`](https://KmBioChemo.github.io/RNAflow/reference/fig_eigengene.md)
  : Module eigengene profile
- [`fig_module_network()`](https://KmBioChemo.github.io/RNAflow/reference/fig_module_network.md)
  : Module co-expression network
- [`enrich_modules()`](https://KmBioChemo.github.io/RNAflow/reference/enrich_modules.md)
  : Enrich every module against a pathway database
- [`fig_module_enrichment()`](https://KmBioChemo.github.io/RNAflow/reference/fig_module_enrichment.md)
  : Module x pathway enrichment dotplot

## Activity inference

Transcription-factor (CollecTRI) and pathway (PROGENy) activity via
decoupleR.

- [`get_tf_network()`](https://KmBioChemo.github.io/RNAflow/reference/get_tf_network.md)
  : Fetch a transcription-factor regulon network (CollecTRI)
- [`get_pathway_network()`](https://KmBioChemo.github.io/RNAflow/reference/get_pathway_network.md)
  : Fetch a pathway-responsive-gene network (PROGENy)
- [`run_activity()`](https://KmBioChemo.github.io/RNAflow/reference/run_activity.md)
  : Infer regulator / pathway activity from a DE contrast
- [`fig_activity_bar()`](https://KmBioChemo.github.io/RNAflow/reference/fig_activity_bar.md)
  : Diverging bar chart of top activity scores

## AI interpretation

Build a prompt from a contrast and ask the Claude API to interpret it.

- [`build_interpret_prompt()`](https://KmBioChemo.github.io/RNAflow/reference/build_interpret_prompt.md)
  : Build the interpretation prompt from a contrast (+ optional
  enrichment)
- [`call_claude()`](https://KmBioChemo.github.io/RNAflow/reference/call_claude.md)
  : Call the Anthropic Claude Messages API
- [`interpret_results()`](https://KmBioChemo.github.io/RNAflow/reference/interpret_results.md)
  : Interpret a DE contrast (+ enrichment) with Claude

## Projects & reproducibility

Save/load sessions, export a Methods script and an HTML report.

- [`save_project()`](https://KmBioChemo.github.io/RNAflow/reference/save_project.md)
  : Save a project to disk
- [`load_project()`](https://KmBioChemo.github.io/RNAflow/reference/load_project.md)
  : Load a project from disk
- [`generate_r_script()`](https://KmBioChemo.github.io/RNAflow/reference/generate_r_script.md)
  : Generate a reproducible R script for an analysis
- [`generate_methods_text()`](https://KmBioChemo.github.io/RNAflow/reference/generate_methods_text.md)
  : Generate a Methods paragraph for an analysis
- [`build_report_html()`](https://KmBioChemo.github.io/RNAflow/reference/build_report_html.md)
  : Build a standalone HTML report for an analysis session

## Shiny modules

UI + server pairs, one per feature.

- [`mod_activity_ui()`](https://KmBioChemo.github.io/RNAflow/reference/mod_activity.md)
  [`mod_activity_server()`](https://KmBioChemo.github.io/RNAflow/reference/mod_activity.md)
  : Activity inference module
- [`mod_ai_ui()`](https://KmBioChemo.github.io/RNAflow/reference/mod_ai.md)
  [`mod_ai_server()`](https://KmBioChemo.github.io/RNAflow/reference/mod_ai.md)
  : AI interpretation module
- [`mod_compare_ui()`](https://KmBioChemo.github.io/RNAflow/reference/mod_compare.md)
  [`mod_compare_server()`](https://KmBioChemo.github.io/RNAflow/reference/mod_compare.md)
  : Multi-contrast comparison module
- [`mod_data_ui()`](https://KmBioChemo.github.io/RNAflow/reference/mod_data.md)
  [`mod_data_server()`](https://KmBioChemo.github.io/RNAflow/reference/mod_data.md)
  : Data input module
- [`mod_de_ui()`](https://KmBioChemo.github.io/RNAflow/reference/mod_de.md)
  [`mod_de_server()`](https://KmBioChemo.github.io/RNAflow/reference/mod_de.md)
  : Differential expression module
- [`mod_enrich_ui()`](https://KmBioChemo.github.io/RNAflow/reference/mod_enrich.md)
  [`mod_enrich_server()`](https://KmBioChemo.github.io/RNAflow/reference/mod_enrich.md)
  : Functional enrichment module
- [`mod_heatmap_ui()`](https://KmBioChemo.github.io/RNAflow/reference/mod_heatmap.md)
  [`mod_heatmap_server()`](https://KmBioChemo.github.io/RNAflow/reference/mod_heatmap.md)
  : Heatmap module
- [`mod_pca_ui()`](https://KmBioChemo.github.io/RNAflow/reference/mod_pca.md)
  [`mod_pca_server()`](https://KmBioChemo.github.io/RNAflow/reference/mod_pca.md)
  : PCA module
- [`mod_project_ui()`](https://KmBioChemo.github.io/RNAflow/reference/mod_project.md)
  [`mod_project_server()`](https://KmBioChemo.github.io/RNAflow/reference/mod_project.md)
  : Project manager module
- [`mod_qc_ui()`](https://KmBioChemo.github.io/RNAflow/reference/mod_qc.md)
  [`mod_qc_server()`](https://KmBioChemo.github.io/RNAflow/reference/mod_qc.md)
  : QC / diagnostics module
- [`mod_report_ui()`](https://KmBioChemo.github.io/RNAflow/reference/mod_report.md)
  [`mod_report_server()`](https://KmBioChemo.github.io/RNAflow/reference/mod_report.md)
  : Reproducibility / report module
- [`mod_volcano_ui()`](https://KmBioChemo.github.io/RNAflow/reference/mod_volcano.md)
  [`mod_volcano_server()`](https://KmBioChemo.github.io/RNAflow/reference/mod_volcano.md)
  : Volcano plot module
- [`mod_wgcna_ui()`](https://KmBioChemo.github.io/RNAflow/reference/mod_wgcna.md)
  [`mod_wgcna_server()`](https://KmBioChemo.github.io/RNAflow/reference/mod_wgcna.md)
  : WGCNA co-expression network module

## Package overview

- [`RNAflow-package`](https://KmBioChemo.github.io/RNAflow/reference/RNAflow-package.md)
  [`RNAflow`](https://KmBioChemo.github.io/RNAflow/reference/RNAflow-package.md)
  : RNAflow: end-to-end bulk RNA-seq analysis platform
- [`%>%`](https://KmBioChemo.github.io/RNAflow/reference/pipe.md) : Pipe
  operator
