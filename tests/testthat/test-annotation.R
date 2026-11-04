test_that("template and validation work", {
  ids <- c("run_a", "run_b", "run_c", "run_d")
  tpl <- annotation_template(ids)
  expect_equal(tpl$sample_id, ids)
  v <- validate_annotation(tpl, ids)
  expect_false(v$valid)                      # groups are empty
  tpl$group <- c("ctrl", "ctrl", "case", "case")
  tpl$sample_name <- c("c1", "c2", "t1", "t2")
  v <- validate_annotation(tpl, ids)
  expect_true(v$valid)
  expect_equal(nrow(v$summary), 2)
  # run mismatch
  v2 <- validate_annotation(tpl, c(ids, "run_e"))
  expect_false(v2$valid)
  expect_true(any(grepl("no annotation row", v2$errors)))
  # duplicated id
  tpl2 <- tpl; tpl2$sample_id[2] <- "run_a"
  expect_false(validate_annotation(tpl2, ids)$valid)
})

test_that("aliases and exclude are normalised", {
  df <- data.frame(Run = c("a", "b"), Condition = c("x", "y"), animal_id = c("1", "2"), exclude = c("yes", "0"))
  n <- labproteomics:::normalise_annotation(df)
  expect_true(all(c("sample_id", "group", "subject_id", "sample_name") %in% names(n)))
  expect_equal(n$exclude, c(TRUE, FALSE))
})

test_that("complete pairs are found", {
  s <- data.frame(sample_id = 1:5, group = c("ctrl", "case", "ctrl", "case", "ctrl"),
                  subject_id = c("A", "A", "B", "B", "C"), stringsAsFactors = FALSE)
  expect_equal(sort(complete_pairs(s, "case", "ctrl")), c("A", "B"))
  expect_equal(complete_pairs(s, "case", "ctrl", pair_by = "missing"), character())
})

test_that("contrast_spec handles pairing", {
  expect_equal(contrast_spec("a", "b")$pairing, "unpaired")
  expect_equal(contrast_spec("a", "b", "subject_id")$pairing, "paired")
  expect_null(contrast_spec("a", "b", "none")$pair_by)
  expect_error(contrast_spec("a", "a"))
})
