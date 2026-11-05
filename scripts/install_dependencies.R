# One-time installation of everything the Lab Proteomics Workbench needs (macOS or Windows).
# Run inside R or RStudio:   source("scripts/install_dependencies.R")
# Then install the workbench itself (see README.md, section "Install").
# No compiler or Rtools is needed: every dependency has a binary on CRAN / Bioconductor, and
# msdap (installed from GitHub) is pure R.
options(repos = c(CRAN = "https://cloud.r-project.org"))
missing <- function(p) p[!vapply(p, requireNamespace, logical(1), quietly = TRUE)]

cran <- c("remotes", "BiocManager", "shiny", "bslib", "DT", "dplyr", "tidyr", "tibble", "readr", "ggplot2", "ggrepel",
          "patchwork", "jsonlite", "rlang", "digest", "msigdbr", "rmarkdown", "knitr", "htmltools", "arrow", "ragg",
          "systemfonts", "readxl", "chromote", "testthat")
if (length(missing(cran))) install.packages(missing(cran))

# gene sets of msigdbr (version 10 and newer keep the full collections in the msigdbdf package)
if (!requireNamespace("msigdbdf", quietly = TRUE)) {
  tryCatch(install.packages("msigdbdf", repos = c("https://igordot.r-universe.dev", "https://cloud.r-project.org")),
           error = function(e) message("msigdbdf could not be installed; msigdbr may offer fewer collections. ", conditionMessage(e)))
}

bioc <- c("limma", "DEqMS", "fgsea", "MSnbase", "ProtGenerics", "vsn", "rmzqc")
if (length(missing(bioc))) BiocManager::install(missing(bioc), update = FALSE, ask = FALSE)

# MS-DAP (peptide-level import, filtering, normalisation, roll-up, DEqMS); from GitHub
if (!requireNamespace("msdap", quietly = TRUE)) remotes::install_github("ftwkoopmans/msdap", upgrade = "never")

still <- missing(c(cran, bioc, "msdap"))
if (length(still)) {
  message("Still missing: ", paste(still, collapse = ", "), ". See https://github.com/ftwkoopmans/msdap for the MS-DAP notes.")
} else {
  message("All dependencies are installed. Now install the workbench package (README, section Install).")
}
