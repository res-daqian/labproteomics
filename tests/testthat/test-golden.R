# Golden tests: the engine must reproduce the locked PhD results on the demo data.
# The demo folder holds the Long-term (stage L) runs only (36 runs; scripts/make_demo_subset.R).
# The Early contrast is checked only when the full 76-run data are present.
demo <- normalizePath(file.path(testthat::test_path(), "..", "..", "demo_data"), mustWork = FALSE)
report <- file.path(demo, "report.parquet")

test_that("the engine reproduces the locked PhD contrasts and rescue counts", {
  skip_if_not(file.exists(report), "demo data not present")
  skip_on_cran()
  ds <- import_quant(report)
  ann <- read_annotation(file.path(demo, "demo_sample_annotation.csv"))
  ds <- attach_annotation(ds, ann)
  expect_true(nrow(ds$samples) %in% c(36, 76))

  if (all(c("LIMA_E", "LIM_E") %in% ds$samples$group)) {
    r1 <- run_contrast(ds, contrast_spec("LIMA_E", "LIM_E", pair_by = "subject_id"), list(limma_check = FALSE))
    expect_equal(r1$n_tested, 5963)
    expect_equal(r1$n_hit_p, 36)
    expect_equal(r1$n_hit_q, 0)
    expect_equal(r1$n_pairs, 10)
  }

  excl <- ds$samples$sample_id[ds$samples$subject_id %in% c("L_CTRL_4", "L_CTRL_BLK")]
  expect_length(excl, 4)
  r2 <- run_contrast(ds, contrast_spec("LIM_L", "CNT_L"), exclude_samples = excl)
  expect_equal(c(r2$n_tested, r2$n_hit_p, r2$n_hit_q), c(6044, 1985, 1915))
  r3 <- run_contrast(ds, contrast_spec("LIMA_L", "LIM_L", pair_by = "subject_id"), exclude_samples = excl)
  expect_equal(c(r3$n_tested, r3$n_hit_p, r3$n_hit_q), c(6193, 47, 0))
  r4 <- run_contrast(ds, contrast_spec("CNTA_L", "CNT_L", pair_by = "subject_id"), list(limma_check = FALSE), exclude_samples = excl)
  expect_equal(c(r4$n_tested, r4$n_hit_p, r4$n_pairs), c(6223, 303, 6))

  res <- rescue_classify(r2, r3, r4)
  expect_equal(res$counts$n[res$counts$tier == "Tier 1"], 15)
  expect_equal(res$counts$n[res$counts$tier == "Tier 2"], 180)
  expect_equal(sum(res$table$drug_only_flag & res$table$tier %in% c("Tier 1", "Tier 2")), 2)

  # Additions of 2026-11-05 (not thesis numbers; regression values on the demo data).
  # The leave-one-out refit on all runs must reproduce the MS-DAP DEqMS p-values exactly.
  expect_lt(r3$loo$full_check, 1e-8)
  expect_lt(r2$loo$full_check, 1e-8)
  expect_equal(r3$loo$n_fragile + r3$loo$n_robust, 47)
  expect_equal(r3$loo$unit_type, "subject")
  expect_equal(r2$loo$unit_type, "run")
  expect_equal(r2$detection$n_candidates, 83)
  expect_true(all(!r2$detection$table$in_main_table))
  expect_equal(r2$n_hit_treat, 687)
  expect_equal(r3$n_hit_treat, 0)
})
