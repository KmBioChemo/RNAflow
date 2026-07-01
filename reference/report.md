# Standalone HTML analysis report

Builds a self-contained HTML report of an RNAflow session using
htmltools – figures are embedded as base64 data URIs, so the file needs
no external assets and no pandoc / Quarto toolchain. Fast, dependable
figures (DE volcanoes, cross-contrast comparison) are rendered inline;
the heavier enrichment / network steps are documented via the
reproducible script (see
[`generate_r_script()`](https://KmBioChemo.github.io/RNAflow/reference/generate_r_script.md)).
