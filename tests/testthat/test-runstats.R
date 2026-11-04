fake_stats <- function(n = 20) {
  set.seed(2)
  tibble::tibble(
    sample_id = sprintf("run_%02d", seq_len(n)), file_name = sprintf("C:\\data\\run_%02d.wiff", seq_len(n)),
    Precursors.Identified = round(stats::rnorm(n, 60000, 800)), Proteins.Identified = round(stats::rnorm(n, 6900, 40)),
    FWHM.RT = stats::rnorm(n, 0.075, 0.001), Median.Mass.Acc.MS2 = stats::rnorm(n, 3.7, 0.05),
    Average.Missed.Tryptic.Cleavages = stats::rnorm(n, 0.1, 0.003)
  )
}

test_that("run scorecard flags a weak run and ignores unknown columns", {
  st <- fake_stats()
  st$Precursors.Identified[5] <- 30000
  st$Average.Missed.Tryptic.Cleavages[7] <- 0.4
  sc <- run_scorecard(st)
  expect_equal(sc$params$strata, "all runs")
  expect_true(sc$table$n_flags[sc$table$sample_id == "run_05"] >= 1)
  expect_true(grepl("Precursors", sc$table$flagged_metrics[sc$table$sample_id == "run_05"]))
  expect_true(grepl("Missed cleavages", sc$table$flagged_metrics[sc$table$sample_id == "run_07"]))
  expect_true(all(c("rules", "systematic", "long") %in% names(sc)))
  expect_s3_class(scorecard_plot(sc), "ggplot")
})

test_that("flags are computed within batch when batches are large enough", {
  st <- fake_stats()
  samples <- data.frame(sample_id = st$sample_id, shortname = st$sample_id, group = rep(c("a", "b"), 10),
                        batch = rep(c("x", "y"), each = 10), run_order = 1:20, stringsAsFactors = FALSE)
  st$Precursors.Identified[1:10] <- st$Precursors.Identified[1:10] - 20000   # batch x is systematically lower
  sc <- run_scorecard(st, samples)
  expect_equal(sc$params$strata, "batch")
  expect_equal(sum(sc$table$flag_Precursors.Identified), 0)                  # no single run leaves its batch
  expect_true("Precursors identified" %in% sc$summary$batch_differences)
})

test_that("stats and log files are found next to a report", {
  d <- tempfile(); dir.create(d)
  writeLines("x", file.path(d, "report.parquet")); writeLines(c("File.Name\tPrecursors.Identified", "a.raw\t1"), file.path(d, "report.stats.tsv"))
  expect_equal(basename(find_diann_stats(file.path(d, "report.parquet"))), "report.stats.tsv")
  expect_true(is.na(find_diann_log(file.path(d, "report.parquet"))))
  st <- read_diann_stats(file.path(d, "report.stats.tsv"))
  expect_equal(st$sample_id, "a")
})
