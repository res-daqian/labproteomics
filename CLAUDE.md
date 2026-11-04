# CLAUDE.md — gui_development_proteomics_workflow

Stable rules for the agent. Moving state (what works, what is next, demo status) lives in
`PROGRESS.md` — read it first each session.

## What this project is

A **locally run Shiny application** (R package `labproteomics`, working name "Lab Proteomics
Workbench") that packages the PhD label-free proteomics workflow for lab colleagues who do not
code. Requested by the supervisor after the 2026-11-04 meeting; first demo at the lab meeting on
Friday 2026-11-06. The analysis methods are those of the PhD project
(`../phd-thesis-atropine`, locked pipeline at `data_locked/pipeline/`, commit c21b4b1):
MS-DAP import → filter → vwmb normalisation → MaxLFQ roll-up → DEqMS, limma-trend concordance,
Hallmark GSEA on ortholog-mapped ranks, and the two-contrast rescue classification (Tier 1 / Tier 2).

## Isolation from the PhD project (never break)

- Never write into `../phd-thesis-atropine`. Demo data were cloned read-only into `demo_data/`
  (gitignored). If a PhD-side export is needed, create it here.
- Functions were **copied and refactored** from the PhD scripts, never sourced across repositories.
  `PROVENANCE.md` maps every engine function to its source script and commit. Keep it current.
- The PhD results stay the reference: the golden tests in `tests/testthat/test-golden.R` must
  reproduce the locked PhD numbers on the demo data (36 Long-term runs; LIM_L vs CNT_L with the two
  excluded animals: 6044 / 1985 / 1915; LIMA_L vs LIM_L: 6193 / 47; CNTA_L vs CNT_L: 6223 / 303;
  Long-term rescue: Tier 1 = 15, Tier 2 = 180). With the full 76-run data the Early contrast
  (LIMA_E vs LIM_E: 5963 / 36) is checked as well. Run them before any commit that touches
  `R/import.R`, `R/qc.R`, `R/dea.R` or `R/rescue.R`.
- Demo data are rebuilt only with `scripts/make_demo_subset.R` followed by
  `scripts/make_demo_annotation.R`. The subset keeps the original row order and only the peptides
  MS-DAP keeps on the full import; both are needed to reproduce the thesis numbers (see
  `PROVENANCE.md`, section "Demo data").

## Architecture rules

- **Engine first, interface second.** Everything the interface does is a call to a plain R function
  in `R/` (`import_quant`, `attach_annotation`, `qc_global`, `run_contrast`, `rescue_classify`,
  `ranked_genes` + `gsea_hallmark`, `write_run_record`). No statistics in Shiny code.
- **Every analysis choice is a parameter** shown in the interface and written to the run record
  (`parameters.json`). No hidden defaults. The PhD defaults are defaults, not constants.
- **Every run writes a run record** (`runs/<timestamp>_<label>/`): input manifest with MD5,
  annotation, parameters, package versions, session info, tables, figures (PDF + PNG), HTML report.
- **Honesty in the interface:** a hit is always stated with its verbatim gate
  (`|log2FC| >= log2(1.2) (= 0.263) AND p < 0.05`); Tier 1 / Tier 2 definitions are printed
  next to the counts; the drug-only contrast only annotates, never filters; unpaired contrasts
  are labelled as such.
- **Figures:** black-grey-red house style (`theme_lab_proteomics`, `group_colours`), thin axes, no grid.
- **Language:** interface and documents in plain ESL-safe English (no idioms, one idea per
  sentence, define abbreviations). Chat with Daqian in Chinese. `PROGRESS.md` is Chinese.
- **Cross-platform:** Mac and Windows. No shell-only tricks in R code; paths via `file.path`;
  no dependency that lacks a Windows binary.
- **Synchronous computation with `withProgress` for now** (single user, each step under one
  minute on the demo data). Async (`ExtendedTask` + `mirai`) is a planned upgrade, not a rewrite.

## Autonomy

- Act then report: commits at coherent checkpoints (2–4 plain sentences for Daqian, a git
  beginner), package installs into the user library (document what and why), mechanics.
- Discuss first: new analysis methods, changes to the rescue definitions or default gates,
  anything shown to the lab or the supervisor.
- Opt-in only: `git push`, destructive git operations.

## Where things are

- `R/` engine (`utils`, `theme`, `annotation`, `import`, `qc`, `dea`, `rescue`, `enrichment`,
  `report`) and interface (`app.R`, `mod_*.R`). `inst/report/run_report.Rmd` report template.
- `tests/testthat/` unit tests + golden tests. `demo_data/` demo inputs (gitignored large files).
- `runs/` run records (gitignored). `scripts/` helpers (`launch_app.R`, `make_demo_annotation.R`).
- `docs/` demo script and notes. `README.md` user guide (Mac / Windows install and start).
