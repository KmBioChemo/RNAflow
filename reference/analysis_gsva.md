# Per-sample gene-set scoring (GSVA / ssGSEA)

Unlike GSEA/ORA (which score a whole contrast), GSVA and ssGSEA assign
every sample its own enrichment score for each gene set, turning a genes
x samples matrix into a sets x samples matrix. That per-sample signature
matrix can be clustered, correlated with phenotype, or fed to downstream
models. Wraps GSVA (a heavy Bioconductor dependency, guarded).
