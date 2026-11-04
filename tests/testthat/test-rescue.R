fake_result <- function(label, log2FC, pvalue) {
  n <- length(log2FC)
  list(spec = list(label = label),
       table = tibble::tibble(protein_id = paste0("P", seq_len(n)), gene_symbol = paste0("G", seq_len(n)),
                              description = "", log2FC = log2FC, pvalue = pvalue, qvalue = pvalue * 2))
}

test_that("rescue tiers follow the locked definitions", {
  disease   <- fake_result("dz", log2FC = c( 1.0,  1.0,  1.0, 0.1,  1.0), pvalue = c(0.01, 0.01, 0.01, 0.5, 0.01))
  treatment <- fake_result("tx", log2FC = c(-0.9, -0.9, -0.2, 0.0, -0.9), pvalue = c(0.01, 0.20, 0.30, 0.9, 0.04))
  res <- rescue_classify(disease, treatment)
  tiers <- stats::setNames(res$table$tier, res$table$protein_id)
  expect_equal(unname(tiers["P1"]), "Tier 1")        # both significant, opposite
  expect_equal(unname(tiers["P2"]), "Tier 2")        # treatment NS, residual 0.1 < 0.263
  expect_equal(unname(tiers["P3"]), "Disease-only")  # residual 0.8 too large
  expect_equal(unname(tiers["P4"]), "Not perturbed")
  expect_equal(unname(tiers["P5"]), "Tier 1")
  expect_equal(res$counts$n[res$counts$tier == "Tier 1"], 2)
  expect_match(res$gate, "log2\\(1.2\\)")
})

test_that("drug-only contrast only annotates", {
  disease   <- fake_result("dz", log2FC = c(1.0, 1.0), pvalue = c(0.01, 0.01))
  treatment <- fake_result("tx", log2FC = c(-0.9, -0.9), pvalue = c(0.01, 0.01))
  drug      <- fake_result("drug", log2FC = c(-0.9, 0.0), pvalue = c(0.01, 0.9))
  res <- rescue_classify(disease, treatment, drug)
  expect_equal(res$table$tier, c("Tier 1", "Tier 1"))
  expect_equal(res$table$drug_only_flag, c(TRUE, FALSE))
})
