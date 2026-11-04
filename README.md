# Lab Proteomics Workbench (working name)

A locally run Shiny application for label-free quantitative proteomics. It packages the
analysis chain of the retina atropine PhD project so that colleagues can run the same
analysis without writing code:

1. **Data** — import DIA-NN (`report.parquet` / `report.tsv`), Spectronaut or FragPipe output;
   attach a sample annotation table (groups, subjects for paired designs, batches).
2. **Quality control** — MS-DAP filtering, vwmb normalisation and MaxLFQ roll-up; missing values,
   intensity centring, PCA, sample correlation, acquisition drift, batch structure,
   within-subject consistency, MAD-based outlier flags; exclude runs.
3. **Differential analysis** — build contrasts (paired or unpaired), set the gate
   (fold change, p or q), run DEqMS with a limma-trend concordance check; volcano plot,
   per-protein abundance plot, result table, CSV download.
4. **Rescue** — classify proteins by whether a treatment reverses a disease effect
   (Tier 1 / Tier 2 with the definitions printed next to the counts).
5. **Enrichment** — Hallmark gene-set enrichment (fgsea) on the model t-statistic, with an
   optional ortholog table (for example guinea pig to human).
6. **Run record** — every analysis writes a folder with the input manifest (MD5), annotation,
   all parameters, package versions, session information, tables, figures (PDF and PNG) and an
   HTML report.

The statistical engine is a set of plain R functions in `R/`; the interface only calls them.
`PROVENANCE.md` lists the origin of every function.

## Requirements

- R 4.4 or newer (R 4.5 tested) on macOS or Windows.
- RStudio is recommended (it provides pandoc for the HTML report).
- R packages: `msdap` (from GitHub, with its Bioconductor dependencies), `shiny`, `bslib`, `DT`,
  `dplyr`, `tidyr`, `readr`, `ggplot2`, `ggrepel`, `patchwork`, `limma`, `DEqMS`, `fgsea`,
  `msigdbr`, `rmarkdown`, `digest`, `jsonlite`, `rlang`. Optional: `readxl` (Excel annotation).

## Install (first time)

```r
install.packages(c("devtools", "BiocManager"))
BiocManager::install(c("limma", "DEqMS", "fgsea", "MSnbase", "ProtGenerics", "vsn"))
devtools::install_github("ftwkoopmans/msdap")        # follow the MS-DAP page if a dependency fails
install.packages(c("shiny", "bslib", "DT", "dplyr", "tidyr", "readr", "ggplot2", "ggrepel",
                   "patchwork", "msigdbr", "rmarkdown", "digest", "jsonlite", "rlang", "readxl"))
devtools::install("path/to/gui_development_proteomics_workflow")   # this folder
```

MS-DAP installation notes for Windows and macOS: https://github.com/ftwkoopmans/msdap

## Start the application

- macOS: double-click `launch_mac.command` (first time: right-click, Open).
- Windows: double-click `launch_windows.bat`.
- From R or RStudio:

```r
labproteomics::run_app(demo_dir = "path/to/demo_data", runs_dir = "path/to/runs")
```

The application opens in the default browser. Close the R process to stop it.

## Sample annotation

| column | required | meaning |
|---|---|---|
| `sample_id` | yes | run identifier exactly as in the search-engine output (file name without extension); the Data page offers a template pre-filled with these names |
| `sample_name` | yes | short label for figures and tables |
| `group` | yes | experimental group; every contrast compares two groups |
| `subject_id` | no | animal or subject; needed for paired designs (for example the two eyes of one animal) |
| `batch` | no | processing batch (for example digestion date) |
| `run_order` | no | acquisition order, used for the drift check |
| `exclude` | no | TRUE leaves the run out of every analysis |

Extra columns are kept and written to the run record.

## Defaults

The defaults are the locked settings of the retina atropine study: at least two peptides per
protein, detection in at least half of the runs of a group, vwmb normalisation, MaxLFQ roll-up,
DEqMS, gate `|log2FC| >= log2(1.2) (= 0.263) AND p < 0.05`. Every setting can be changed in the
interface and is written to the run record.

## Status (2026-11-04)

First version for the lab-meeting demo. DIA-NN input is tested with real data: the full 76-run report
of the retina atropine study during development, and the 36 Long-term runs that form the demo data
(the demo reproduces the thesis numbers). Spectronaut and FragPipe imports call the MS-DAP importers
but have not been tested with data. Computation runs in the foreground with a progress bar; the
global QC takes about half a minute for 36 runs. Enrichment covers the Hallmark collection only.
The demo data are not part of the package or the repository; `scripts/make_demo_subset.R` followed
by `scripts/make_demo_annotation.R` rebuilds them from the PhD inputs (see `PROVENANCE.md`).
