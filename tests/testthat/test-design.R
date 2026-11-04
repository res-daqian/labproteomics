make_samples <- function() {
  data.frame(
    sample_id = sprintf("run_%02d", 1:12), sample_name = sprintf("s%02d", 1:12),
    group = rep(c("ctrl", "case", "treated"), each = 4),
    subject_id = c(paste0("A", 1:4), paste0("A", 1:4), paste0("B", 1:4)),
    batch = c(rep("b1", 4), rep("b2", 4), rep("b2", 4)),
    run_order = c(1, 3, 5, 7, 2, 4, 6, 8, 9, 10, 11, 12),
    stringsAsFactors = FALSE
  )
}

test_that("design check finds confounded and crossed group pairs", {
  dc <- design_check(make_samples())
  gp <- dc$group_pairs
  pair <- function(a, b) gp[(gp$group_a == a & gp$group_b == b) | (gp$group_a == b & gp$group_b == a), ]
  expect_equal(pair("ctrl", "case")$batch_status, "confounded")
  expect_equal(pair("case", "treated")$batch_status, "same single batch")
  expect_equal(pair("ctrl", "case")$complete_pairs, 4)
  m <- design_messages(dc)
  expect_true(any(m$level == "danger"))
  expect_true(any(grepl("No pooled QC", m$text)))
  expect_true(dc$has$batch && dc$has$subject && dc$has$run_order)
})

test_that("design check recognises QC injections and works without optional columns", {
  s <- make_samples()
  s$group[12] <- "QC_pool"
  dc <- design_check(s)
  expect_equal(nrow(dc$qc_runs), 1)
  expect_true(any(grepl("pooled QC or standard injection", design_messages(dc)$text)))
  minimal <- data.frame(sample_id = letters[1:4], group = c("a", "a", "b", "b"), stringsAsFactors = FALSE)
  dc2 <- design_check(minimal)
  expect_null(dc2$batch_table)
  expect_true(any(grepl("No batch column", design_messages(dc2)$text)))
})
