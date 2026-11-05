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

- Never write into `../phd-thesis-atropine`. If a PhD-side export is needed, create it here.
- **No laboratory data ship with the workbench.** The demo (`demo_data/`) is the public PRIDE data
  set PXD079891 (HUVEC, MANCR knockdown vs control, 8 donors, 16 runs, DIA-NN 1.9, CC0), fetched by
  `scripts/get_demo_data.R`; the parquet and the optional FASTA are gitignored, the small files
  (stats, log, SDRF, annotation, settings, README) are committed. Daqian decided this on
  2026-11-05 ("放弃用我自己的数据").
- The thesis subset (36 Long-term runs, Daqian's own data) lives in `golden_data/` (fully
  gitignored, read-only copies). It is used by the golden tests and, on Daqian's Mac only, as the
  data of the Friday demo (`launch_mac_thesis_data.command` = `Rscript scripts/launch_app.R
  golden_data`; Daqian decided on 2026-11-05 evening to present his own data). It is rebuilt with
  `scripts/make_golden_subset.R` followed by `scripts/make_golden_annotation.R`. A folder used by
  the demo button needs `demo_settings.json` (label, tissue, species); the one in `golden_data/`
  is local too.
- Functions were **copied and refactored** from the PhD scripts, never sourced across repositories.
  `PROVENANCE.md` maps every engine function to its source script and commit. Keep it current.
- The PhD results stay the reference: the golden tests in `tests/testthat/test-golden.R` must
  reproduce the locked PhD numbers on `golden_data/` (36 Long-term runs; LIM_L vs CNT_L with the two
  excluded animals: 6044 / 1985 / 1915; LIMA_L vs LIM_L: 6193 / 47; CNTA_L vs CNT_L: 6223 / 303;
  Long-term rescue: Tier 1 = 15, Tier 2 = 180). With the full 76-run data the Early contrast
  (LIMA_E vs LIM_E: 5963 / 36) is checked as well. Run them with `NOT_CRAN=true` before any commit
  that touches `R/import.R`, `R/qc.R`, `R/dea.R` or `R/rescue.R`.
- The golden subset keeps the original row order and only the peptides MS-DAP keeps on the full
  import; both are needed to reproduce the thesis numbers (see `PROVENANCE.md`, "Golden data").
- The git tag `demo-2026-11-06-safe` marks the first demo-ready version (2026-11-04). Keep it.

## Architecture rules

- **Engine first, interface second.** Everything the interface does is a call to a plain R function
  in `R/` (`import_quant`, `attach_annotation`, `design_check`, `run_scorecard`, `marker_scores`,
  `qc_global`, `run_contrast` (with `detection_differences` and `loo_stability`), `rescue_classify`,
  `ranked_genes` + `gsea_collection` / `ora_collection`, `export_adapters`, `write_run_record`).
  No statistics in Shiny code.
- **Every analysis choice is a parameter** shown in the interface and written to the run record
  (`parameters.json`). No hidden defaults. The PhD defaults are defaults, not constants.
- **Every run writes a run record** (`runs/<timestamp>_<label>/`): input manifest with MD5,
  annotation, parameters, package versions, session info, tables, figures (PDF + PNG), HTML report.
- **Honesty in the interface:** a hit is always stated with its verbatim gate
  (`|log2FC| >= log2(1.2) (= 0.263) AND p < 0.05`, or the q / treat variant); Tier 1 / Tier 2
  definitions are printed next to the counts; the drug-only contrast only annotates, never filters;
  unpaired contrasts are labelled as such; leave-one-out fragility, detection-only differences and
  treat columns describe the hits and never change them; the ORA background is always the
  quantified proteins of the contrast.
- **Figures:** Daqian's palettes (2026-11-05; constants in `R/theme.R`: qualitative
  `LP_QUALITATIVE`, diverging `LP_DIVERGING`, sequential `LP_SEQ_COOL` / `LP_SEQ_WARM`), Arial,
  left and bottom hairline axes only, no grid, no panel box, narrow bars, open circles for samples.
  Red = case / increase / flag, teal = control / decrease / main series. PDF via cairo, PNG via ragg.
- **Interface:** teal `#015069` as the working colour, red `#c71e1d` only for alerts and hits;
  body font Avenir Next / Segoe UI (system stack), headings and tile values in Charter / Palatino /
  Georgia (system serif stack); no web fonts are downloaded. The Rescue page is shown only when
  the attached design has three or more groups (`nav_show` / `nav_hide` in `app_server`).
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
- Opt-in only: destructive git operations (force push, history rewrite, branch deletion). Pushing
  `main` after each commit is routine since 2026-11-05 (Daqian's instruction).

## Where things are

- `R/` engine (`utils`, `theme`, `annotation`, `import`, `design`, `runstats`, `markers`, `qc`,
  `dea`, `dea_extras`, `rescue`, `enrichment`, `export`, `report`, `qc_report`) and interface
  (`app.R` with the theme and the building blocks `stat_tile`, `page_header`, `msg_list`;
  `mod_*.R`). `inst/report/run_report.Rmd` run report; `inst/report/qc_report.Rmd` + `.css`
  stand-alone QC report (HTML, printed to PDF by chromote + Chrome; mzQC via rmzqc, local CV only,
  single `MzQCcvParameter` for the batch label because a list nests and fails the schema).
- `tests/testthat/` unit tests + golden tests. `demo_data/` public demo (large files gitignored;
  `demo_data/README.md`). `golden_data/` thesis subset for the golden tests (gitignored).
- `runs/` run records (gitignored). `scripts/` helpers (`launch_app.R`, `get_demo_data.R`,
  `make_golden_subset.R`, `make_golden_annotation.R`).
- `docs/` demo scripts (`demo_script_CN.md` = thesis data, `demo_script_public_data_CN.md` =
  public data), QC primer (`qc_primer_CN.md`), parameters and stack explained
  (`parameters_and_stack_CN.md`, private), lab-meeting handout for colleagues
  (`lab_meeting_handout.md`, English) and talk outline (`lab_meeting_talk_CN.md`). `README.md`
  user guide; `scripts/install_dependencies.R` one-shot dependency installer for colleagues.
- Public repository since 2026-11-05 evening: https://github.com/res-daqian/labproteomics
  (remote `origin`, branch `main`, tag `demo-2026-11-06-safe` pushed; Daqian created it himself
  with `gh repo create` because the harness blocks publication actions). Since 2026-11-05 evening
  Daqian wants every commit pushed right away ("以后每次提交都顺手推送"): commit, then
  `git push origin main`. Destructive git operations (history rewrite, force push) stay opt-in.
