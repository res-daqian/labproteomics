# PROVENANCE — where every engine function comes from

Source repository: `../phd-thesis-atropine`, commit `c21b4b1` (2026-11-02), folder `data_locked/pipeline/`.
Functions were copied and refactored; nothing is sourced across repositories.

| labproteomics function (file) | PhD source | Change made here |
|---|---|---|
| `%||%`, `first_token`, `ensure_dir`, `write_json_pretty`, `coalesce_char_columns`, `fc_log2`, `gate_label` (`R/utils.R`) | `common.R` | vectorised `first_token`; `fc_log2(fc)` takes the fold change as an argument; `gate_label` knows the p, q and treat gates |
| `theme_lab_proteomics`, palette constants, `group_colours` (`R/theme.R`) | `scripts/theme_thesis.R` | `group_colours` falls back to case = red, control = black, others grey for any group names |
| `read_annotation`, `validate_annotation`, `annotation_template`, `complete_pairs` (`R/annotation.R`) | `common.R` (`load_project_metadata`, `translate_metadata_for_msdap`, `pair_summary_table`) + the retired Python metadata check called by `run_stage_01_qc.R` | generic columns (`sample_id`, `sample_name`, `group`, `subject_id`, `batch`, `run_order`, `exclude`); aliases accepted; checks implemented in R |
| `import_quant`, `attach_annotation`, `subset_samples`, `protein_annotations`, `resolve_intensity_column` (`R/import.R`) | `common.R` (`import_msdap_project`, `subset_dataset_samples`, `standardize_protein_metadata`, `resolve_intensity_column`) | format detection; Spectronaut / FragPipe wrappers (untested); optional FASTA import for descriptions |
| `qc_global`, `qc_plots`, `flag_by_mad*`, `pair_consistency_table` (`R/qc.R`) | `run_stage_01_qc.R` | strata for MAD flags = group (when every group has >= 6 runs) instead of stage x pair_family; pair consistency over `subject_id`; batch / run-order plots optional; black-grey-red |
| `dea_defaults`, `contrast_spec`, `run_contrast`, `limma_concordance`, `dea_summary`, `volcano_plot`, `protein_profile_plot` (`R/dea.R`) | `run_stage_02_dea_enhanced.R` (`run_single_contrast`, `run_limma_contrast`), `common.R` (`add_project_contrast`, `locked_hit`) | parameters exposed; gate on p, q or limma treat; complete pairs enforced; per-run profile plot; quantification counts per group; treat columns from the limma-trend fit; detection-only and leave-one-out results attached |
| `rescue_classify`, `rescue_plot` (`R/rescue.R`) | `run_stage_03b_rescue_v2.R` (2-case framework, 2026-05-12) | identical logic; user-facing tiers only; `Disease-only` and `Not perturbed` labels added; `gate` argument; fragility columns carried from the contrasts |
| `read_ortholog_map`, `mapped_symbols`, `ranked_genes`, `gsea_hallmark` (`R/enrichment.R`) | `run_stage_03a_enrichment_v2.R` (ranking, ortholog join, deduplication) | fgsea on MSigDB; the ranking and mapping are unchanged |
| `new_run_dir`, `write_run_record`, `render_run_report` (`R/report.R`, `inst/report/run_report.Rmd`) | `common.R` (`ensure_run_scaffold`: `00_input_freeze` pattern) | one folder per run with MD5 manifest, parameters, versions, tables, figures, adapters, HTML report |

## New in labproteomics (2026-11-05; no PhD source)

These were written for the workbench after the discussion of 2026-11-04/05 on quality control.
They do not change any thesis number (the golden tests cover the contrasts and the rescue counts).

| function (file) | what it does | notes on the method |
|---|---|---|
| `design_check`, `design_messages` (`R/design.R`) | batch x group contingency; every pair of groups classified as confounded / crossed / same single batch; Kruskal-Wallis on run order by group; subjects with runs in several groups; pooled QC injections recognised by name | motivated by the digestion-batch / lens-status confounding of the retina study; the PhD stage-01 QC reported the same confounding with PERMANOVA |
| `find_diann_stats`, `find_diann_log`, `read_diann_stats`, `read_diann_log`, `run_scorecard`, `scorecard_plot`, `scorecard_metrics` (`R/runstats.R`) | per-run instrument and identification metrics from DIA-NN `report.stats.tsv`; flags within batch (median +/- k MAD; identifications below 70% of the median; missed cleavages above 0.25); Kruskal-Wallis per metric by batch and by group ("systematic differences") | DIA-NN metric definitions follow the DIA-NN documentation; thresholds are conventions, shown in the interface |
| `marker_panels`, `marker_roles`, `read_marker_panel`, `marker_scores`, `marker_plot` (`R/markers.R`) | share of the total linear ion signal per marker set and run (unnormalised peptide intensities); expected sets flagged when low, neighbour / contamination sets when high (k MAD and a 0.1% floor) | marker lists compiled from tissue-expression knowledge (Human Protein Atlas tissue-enriched genes, standard contaminant lists); to be reviewed per tissue with Daqian |
| `detection_differences` (`R/dea_extras.R`) | proteins detected in >= 80% of one group and <= 20% of the other; Fisher exact test on detection counts (ignores pairing); unnormalised log2 signal | separate from the fold-change table because such proteins have no trustworthy fold change |
| `deqms_fit`, `treat_fit`, `loo_stability` (`R/dea_extras.R`) | DEqMS refit on the fixed protein matrix (mirrors `msdap:::de_deqms`: eBayes without trend, sigma guard, `spectraCounteBayes`) leaving out one subject or run at a time; a hit lost in any refit is fragile | the refit on all runs reproduces the MS-DAP p-values to 1e-13 (checked in the golden tests); the quantification is not repeated |
| `gene_set_collections`, `load_gene_sets`, `gsea_collection`, `ora_collection`, `enrichment_plot`, `ora_plot` (`R/enrichment.R`) | Hallmark, GO BP/CC/MF, Reactome, KEGG MEDICUS, C8 cell types from msigdbr (offline, cached per session); ORA with `fgsea::fora` against the quantified and mapped proteins of the contrast | background choice follows Wijesooriya et al. 2022 (PLoS Comput Biol) |
| `export_adapters` (`R/export.R`) | IPA, iPathwayGuide, STRING and Cytoscape tables per contrast with a manifest | format rules follow the PhD skill `network-and-pathway-adapters` (gene symbol + protein id + log2FC + significance; contrast in every row; rescue tier in node tables) |

Not ported (thesis-specific): stages 03c–03l (profile plots for the PRM panel, validation ranking,
pathway crosstalk, biology themes, evidence matrix, PRM selection).

## Demo data (`demo_data/`; the parquet, the FASTA, the stats file and the log are gitignored)

| file | source (read only) | how it was made |
|---|---|---|
| `report.parquet` (36 runs, 361 MB) | PhD `data_locked/inputs/diann_v231_output_msaccuracy_auto/report.parquet` (DIA-NN 2.3.1, 76 runs) | `scripts/make_demo_subset.R` (2026-11-04): Long-term (stage L) runs only; rows limited to the peptides MS-DAP keeps on the full 76-run import (best charge state per modified sequence, plain sequences mapped to one protein group); original row order preserved. |
| `report.stats.tsv` (36 rows) | PhD `.../report.stats.tsv` (76 rows) | rows of the 36 Long-term runs (2026-11-05) |
| `report.log.txt` | PhD `.../report.log.txt` | copied unchanged (DIA-NN version, command line and settings) |
| `_phd_sample_metadata_source.csv` | PhD `data_locked/inputs/sample_metadata.csv` | same script; stage L rows only |
| `demo_sample_annotation.csv` | derived from the two files above | `scripts/make_demo_annotation.R` |
| `guineapig_rev_and_unrev_2025_12_19.fasta` | PhD `data_locked/inputs/` | copied read-only |
| `guinea_pig_to_human_ortholog_map.csv` | PhD project ortholog map | copied |

Why the subset is built this way (both findings are from 2026-11-04, verified with the golden tests):

1. MS-DAP decides per data set which charge state of a modified sequence is kept and which plain
   sequences are ambiguous (mapped to more than one protein group). A plain 36-run row filter keeps
   peptides the 76-run import had removed, and the Long-term contrasts then differ by a few hits
   (LIM_L vs CNT_L 1962 instead of 1985; CNTA_L vs CNT_L 301 instead of 303; Tier 2 163 instead of 180).
2. MS-DAP's normalisation and roll-up are numerically sensitive to the row order of the input. With
   identical peptides but scrambled rows, protein intensities moved by up to 0.034 log2 and 8
   borderline hits (p between 0.02 and 0.056) changed status in CNTA_L vs CNT_L (309 instead of 303).
   Reading the parquet as one table and subsetting in R keeps the order, and the numbers are exact.
