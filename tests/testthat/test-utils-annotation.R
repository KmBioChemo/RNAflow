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
