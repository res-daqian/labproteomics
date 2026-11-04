fake_ds <- function() {
  genes <- c("RHO", "PDE6B", "CRYAA", "CRYAB", "HBB", "ACTB")
  prot <- tibble::tibble(protein_id = paste0("P", seq_along(genes)), gene_symbols_or_id = genes, fasta_headers = paste("sp|x|", genes, "desc"))
  samples <- sprintf("run_%02d", 1:8)
  set.seed(3)
  pep <- expand.grid(protein_id = prot$protein_id, sample_id = samples, stringsAsFactors = FALSE)
  pep$peptide_id <- paste0(pep$protein_id, "_pep")
  pep$intensity <- stats::rnorm(nrow(pep), 10, 0.2)
  pep$intensity[pep$protein_id %in% c("P3", "P4") & pep$sample_id == "run_04"] <- 16   # lens spike in run 4
  pep$detect <- TRUE
  list(peptides = tibble::as_tibble(pep), proteins = prot,
       samples = tibble::tibble(sample_id = samples, shortname = samples, group = rep(c("a", "b"), 4), run_order = 1:8))
}

test_that("marker shares are computed and a crystallin spike is flagged", {
  ms <- marker_scores(fake_ds(), tissue = "retina")
  expect_true(all(ms$long$share >= 0 & ms$long$share <= 1))
  lens <- ms$long[ms$long$set == "lens", ]
  expect_true(lens$flagged[lens$sample_id == "run_04"])
  expect_equal(sum(lens$flagged), 1)
  expect_true("run_04" %in% ms$summary$runs_flagged)
  expect_equal(ms$panel$n_proteins_in_data[ms$panel$set == "photoreceptor"], 2)
  expect_s3_class(marker_plot(ms), "ggplot")
})

test_that("roles follow the tissue and custom panels are read", {
  r <- marker_roles("tears")
  expect_equal(r$role[r$set == "tear_film"], "expected")
  expect_equal(r$role[r$set == "lens"], "contamination")
  f <- tempfile(fileext = ".csv")
  writeLines(c("set,gene,role", "lens,CRYAA,contamination", "lens,CRYAB,contamination", "house,ACTB,expected"), f)
  p <- read_marker_panel(f)
  expect_equal(nrow(p), 2)
  expect_equal(p$role[p$set == "house"], "expected")
  ms <- marker_scores(fake_ds(), panel = p)
  expect_setequal(unique(ms$long$set), c("lens", "house"))
})
