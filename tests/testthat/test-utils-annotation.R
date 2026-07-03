test_that("organism_info resolves the three supported organisms", {
  expect_equal(organism_info("human")$orgdb, "org.Hs.eg.db")
  expect_equal(organism_info("Mouse")$kegg, "mmu")       # case-insensitive
  expect_equal(organism_info("rat")$reactome, "rat")
  expect_setequal(supported_organisms(), c("human", "mouse", "rat"))
})

test_that("organism_info rejects bad input", {
  expect_error(organism_info("zebrafish"), "Unsupported organism")
  expect_error(organism_info(NULL), "single value")
  expect_error(organism_info(c("human", "mouse")), "single value")
})

test_that("symbols_to_entrez maps known mouse symbols and drops the rest", {
  skip_if_not_installed("org.Mm.eg.db")
  skip_if_not_installed("AnnotationDbi")
  out <- symbols_to_entrez(c("Actb", "Gapdh", "NOT_A_GENE_XYZ"),
                           "mouse", quiet = TRUE)
  expect_true("Actb" %in% names(out))
  expect_false("NOT_A_GENE_XYZ" %in% names(out))
  expect_type(out, "character")
  # ENTREZ IDs are digit strings
  expect_true(all(grepl("^[0-9]+$", out)))
})

test_that("symbols_to_entrez handles empty / all-unmapped input", {
  skip_if_not_installed("org.Mm.eg.db")
  expect_length(symbols_to_entrez(character(0), "mouse", quiet = TRUE), 0)
  expect_length(symbols_to_entrez(c("ZZZ1", "ZZZ2"), "mouse", quiet = TRUE), 0)
})

test_that("guess_id_type distinguishes ensembl / entrez / symbol", {
  expect_equal(guess_id_type(c("ENSG00000141510.3", "ENSG00000012048")), "ensembl")
  expect_equal(guess_id_type(c("ENSMUSG00000059552", "ENSMUSG00000017950")), "ensembl")
  expect_equal(guess_id_type(c("7157", "672", "5290")), "entrez")
  expect_equal(guess_id_type(c("TP53", "BRCA1", "DUSP1")), "symbol")
  expect_equal(guess_id_type(character(0)), "symbol")
})

test_that("map_de_to_symbols converts Ensembl IDs and collapses duplicates", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("AnnotationDbi")
  res <- data.frame(
    gene = c("ENSG00000141510.3", "ENSG00000012048", "ENSG00000141510",
             "ENSG00000000000"),
    baseMean = c(100, 50, 80, 10), log2FoldChange = c(2, -1, 1.5, 0.2),
    stat = c(5, -3, 4, 0.1), pvalue = 1e-3, padj = 1e-2,
    stringsAsFactors = FALSE)
  out <- map_de_to_symbols(res, "human")
  expect_equal(attr(out, "id_converted"), "ensembl")
  expect_true("TP53" %in% out$gene && "BRCA1" %in% out$gene)
  expect_false(anyDuplicated(out$gene) > 0)          # duplicate ENSG collapsed
  expect_equal(out$stat[out$gene == "TP53"], 5)      # kept the more significant row
})

test_that("map_de_to_symbols leaves symbol input unchanged", {
  res <- data.frame(gene = c("TP53", "BRCA1"), stat = c(1, 2),
                    log2FoldChange = c(1, -1), padj = c(0.01, 0.02))
  out <- map_de_to_symbols(res, "human")
  expect_null(attr(out, "id_converted"))
  expect_identical(out$gene, res$gene)
})
