# Lab Proteomics Workbench (working name)

A locally run Shiny application for label-free quantitative proteomics. It packages the
analysis chain of the retina atropine PhD project so that colleagues can run the same
analysis without writing code:

1. **Data** — import DIA-NN (`report.parquet` / `report.tsv`), Spectronaut or FragPipe output;
   attach a sample annotation table (groups, subjects for paired designs, batches, run order);
   design check before any statistics (batch against group, run order, pairing, pooled QC
   injections); DIA-NN version and settings read from `report.log.txt`.
2. **Quality control** — three layers. Run scorecard from the search engine's per-run statistics
   (`report.stats.tsv`: identifications, signal, peak width, mass error, retention-time deviation,
   normalisation instability, missed cleavages; flags within batch, systematic batch differences
   reported separately). Tissue purity and contamination marker panel (photoreceptor, RPE, lens,
   blood, plasma, sclera/choroid ECM, muscle, skin keratins, ocular surface epithelium, tear film,
   Y-linked; roles depend on the tissue; custom panels accepted). Cohort QC after MS-DAP filtering,
   vwmb normalisation and MaxLFQ roll-up: missing values, intensity centring, PCA, sample
   correlation, acquisition drift, batch structure, within-subject consistency, MAD-based outlier
   flags; exclude runs.
3. **Differential analysis** — build contrasts (paired or unpaired), set the gate (fold change
   with p or q, or limma treat against the fold-change threshold), run DEqMS with a limma-trend
   concordance check; detection-only differences (proteins seen in one group only, Fisher exact
   test) listed separately; leave-one-out stability of every hit (hollow points in the volcano
   plot); per-protein abundance plot with quantification counts; result tables, CSV download.
4. **Rescue** — classify proteins by whether a treatment reverses a disease effect
   (Tier 1 / Tier 2 with the definitions printed next to the counts; fragility of each hit carried
   from the leave-one-out refits).
5. **Enrichment** — GSEA (fgsea, ranked by the model t-statistic) or ORA (hypergeometric, with
   the proteins quantified in the contrast as background) on MSigDB Hallmark, GO (BP, CC, MF),
   Reactome, KEGG MEDICUS or cell-type signatures, with an optional ortholog table (for example
   guinea pig to human). Gene sets come from msigdbr and work offline.
6. **Run record** — every analysis writes a folder with the input manifest (MD5), search-engine
   settings, annotation, design check, scorecard, marker panel, all parameters, package versions,
   session information, tables, figures (PDF and PNG), upload tables for IPA, iPathwayGuide,
   STRING and Cytoscape, and an HTML report.

The statistical engine is a set of plain R functions in `R/`; the interface only calls them.
`PROVENANCE.md` lists the origin of every function.

## Requirements

- R 4.4 or newer (R 4.5 tested) on macOS or Windows.
- RStudio is recommended (it provides pandoc for the HTML report).
- R packages: `msdap` (from GitHub, with its Bioconductor dependencies), `shiny`, `bslib`, `DT`,
  `dplyr`, `tidyr`, `readr`, `ggplot2`, `ggrepel`, `patchwork`, `limma`, `DEqMS`, `fgsea`,
  `msigdbr` (version 10 or newer; the gene sets are bundled), `rmarkdown`, `digest`, `jsonlite`,
  `rlang`. Optional: `readxl` (Excel annotation).

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
interface and is written to the run record. Additional outputs of the differential analysis
(detection-only differences with present >= 80% / absent <= 20%, limma treat columns, leave-one-out
stability) do not change the hits; they describe them.

## Quality-control inputs

- The run scorecard needs the search engine's per-run statistics next to the report
  (DIA-NN: `report.stats.tsv`; found automatically). Without it the scorecard card says so and the
  other QC layers still run.
- The search-engine settings are read from `report.log.txt` next to the report when present.
- The marker panel uses the gene symbols of the data (from the report or the FASTA). Choose the
  tissue on the QC page; a custom panel is a CSV with columns `set`, `gene` and optional `role`
  (expected, neighbour, contamination, check).
- Pooled QC or standard injections are recognised by `pool`, `QC`, `HeLa`, `standard` or `blank`
  in the group label or the sample name. Keep them out of the contrasts.

## Status (2026-11-05)

Second version. DIA-NN input is tested with real data: the full 76-run report of the retina
atropine study during development, and the 36 Long-term runs that form the demo data (the demo
reproduces the thesis numbers; golden tests). Spectronaut and FragPipe imports call the MS-DAP
importers but have not been tested with data; the run scorecard reads DIA-NN statistics only.
Computation runs in the foreground with a progress bar; the global QC takes about half a minute
for 36 runs and a contrast 12 to 27 seconds including the leave-one-out refits. The demo data are
not part of the package or the repository; `scripts/make_demo_subset.R` followed by
`scripts/make_demo_annotation.R` rebuilds them from the PhD inputs (see `PROVENANCE.md`).
Not yet done: phenotype association, msqrob2 as a second engine, chick ortholog map, Windows
real-machine test, asynchronous computation.
