test_that("gate labels cover p, q and treat", {
  expect_match(gate_label(1.2, 0.05, "p"), "AND p < 0.05")
  expect_match(gate_label(1.2, 0.05, TRUE), "AND q < 0.05")
  expect_match(gate_label(1.2, 0.05, "treat"), "treat q < 0.05")
  expect_equal(labproteomics:::resolve_gate(list(gate = "", use_adjusted = TRUE)), "q")
  expect_error(labproteomics:::resolve_gate(list(gate = "x")))
})

test_that("detection-only differences are found on a toy dataset", {
  samples <- tibble::tibble(sample_id = sprintf("s%d", 1:10), shortname = sprintf("s%d", 1:10), group = rep(c("ctrl", "case"), each = 5))
  rows <- list()
  add <- function(pid, sids) rows[[length(rows) + 1]] <<- tibble::tibble(protein_id = pid, sample_id = sids, peptide_id = paste0(pid, "_p"), intensity = 10, detect = TRUE)
  add("A", samples$sample_id)                         # everywhere
  add("B", samples$sample_id[6:10])                   # case only
  add("C", samples$sample_id[1:4])                    # control only (4 of 5)
  add("D", samples$sample_id[c(1, 6)])                # sparse
  d <- list(peptides = dplyr::bind_rows(rows), samples = samples, proteins = tibble::tibble(protein_id = c("A", "B", "C", "D"), gene_symbols_or_id = c("A", "B", "C", "D")))
  det <- detection_differences(d, contrast_spec("case", "ctrl"))
  expect_setequal(det$table$protein_id, c("B", "C"))
  expect_equal(det$n_case_only, 1)
  expect_equal(det$n_control_only, 1)
  expect_true(all(det$table$fisher_p <= 0.05))
})

test_that("deqms refit matches a direct DEqMS fit", {
  set.seed(1)
  mat <- matrix(stats::rnorm(200 * 8, 10), nrow = 200, dimnames = list(paste0("P", 1:200), paste0("s", 1:8)))
  mat[1:10, 5:8] <- mat[1:10, 5:8] + 1
  samples <- tibble::tibble(sample_id = colnames(mat), shortname = colnames(mat), group = rep(c("ctrl", "case"), each = 4))
  npep <- stats::setNames(sample(1:10, 200, replace = TRUE), rownames(mat))
  f <- labproteomics:::deqms_fit(mat, samples, contrast_spec("case", "ctrl"), npep)
  design <- stats::model.matrix(~ rep(0:1, each = 4))
  fit <- limma::eBayes(limma::lmFit(mat, design)); fit$count <- npep[rownames(mat)]; fit <- DEqMS::spectraCounteBayes(fit)
  expect_equal(f$pvalue, unname(fit$sca.p[, 2]), tolerance = 1e-12)
  expect_true(all(f$protein_id[1:10] %in% rownames(mat)[1:10]))
})

test_that("rescue uses the requested gate and carries fragility", {
  mk <- function(label, fc, p, fragile) list(spec = list(label = label), table = tibble::tibble(protein_id = c("P1", "P2"), gene_symbol = c("G1", "G2"), description = "",
                                                                                               log2FC = fc, pvalue = p, qvalue = p, treat_q = c(0.01, 0.9), loo_fragile = fragile))
  dz <- mk("dz", c(1, 1), c(0.01, 0.01), c(TRUE, FALSE)); tx <- mk("tx", c(-1, -1), c(0.01, 0.01), c(FALSE, TRUE))
  res <- rescue_classify(dz, tx)
  expect_equal(res$table$disease_fragile, c(TRUE, FALSE))
  expect_equal(res$table$treatment_fragile, c(FALSE, TRUE))
  res_t <- rescue_classify(dz, tx, gate = "treat")
  expect_equal(res_t$table$tier[res_t$table$protein_id == "P1"], "Tier 1")
  expect_equal(res_t$table$tier[res_t$table$protein_id == "P2"], "Not perturbed")
  expect_match(res_t$gate, "treat")
})
