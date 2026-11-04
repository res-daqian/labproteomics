# PROVENANCE — where every engine function comes from

Source repository: `../phd-thesis-atropine`, commit `c21b4b1` (2026-11-02), folder `data_locked/pipeline/`.
Functions were copied and refactored; nothing is sourced across repositories.

| labproteomics function (file) | PhD source | Change made here |
|---|---|---|
| `%||%`, `first_token`, `ensure_dir`, `write_json_pretty`, `coalesce_char_columns`, `fc_log2` (`R/utils.R`) | `common.R` | vectorised `first_token`; `fc_log2(fc)` takes the fold change as an argument |
| `theme_lab_proteomics`, palette constants, `group_colours` (`R/theme.R`) | `scripts/theme_thesis.R` | `group_colours` falls back to case = red, control = black, others grey for any group names |
| `read_annotation`, `validate_annotation`, `annotation_template`, `complete_pairs` (`R/annotation.R`) | `common.R` (`load_project_metadata`, `translate_metadata_for_msdap`, `pair_summary_table`) + the retired Python metadata check called by `run_stage_01_qc.R` | generic columns (`sample_id`, `sample_name`, `group`, `subject_id`, `batch`, `run_order`, `exclude`); aliases accepted; checks implemented in R |
| `import_quant`, `attach_annotation`, `subset_samples`, `protein_annotations`, `resolve_intensity_column` (`R/import.R`) | `common.R` (`import_msdap_project`, `subset_dataset_samples`, `standardize_protein_metadata`, `resolve_intensity_column`) | format detection; Spectronaut / FragPipe wrappers (untested); optional FASTA import for descriptions |
| `qc_global`, `qc_plots`, `flag_by_mad*`, `pair_consistency_table` (`R/qc.R`) | `run_stage_01_qc.R` | strata for MAD flags = group (when every group has >= 6 runs) instead of stage x pair_family; pair consistency over `subject_id`; batch / run-order plots optional; black-grey-red |
| `dea_defaults`, `contrast_spec`, `run_contrast`, `limma_concordance`, `dea_summary`, `volcano_plot`, `protein_profile_plot` (`R/dea.R`) | `run_stage_02_dea_enhanced.R` (`run_single_contrast`, `run_limma_contrast`), `common.R` (`add_project_contrast`, `locked_hit`) | parameters exposed; gate on p or q; complete pairs enforced; per-run profile plot |
| `rescue_classify`, `rescue_plot` (`R/rescue.R`) | `run_stage_03b_rescue_v2.R` (2-case framework, 2026-05-12) | identical logic; user-facing tiers only; `Disease-only` and `Not perturbed` labels added for the remainder |
| `ranked_genes`, `gsea_hallmark`, `gsea_plot`, `read_ortholog_map` (`R/enrichment.R`) | `run_stage_03a_enrichment_v2.R` (ranking, ortholog join, deduplication) | fgsea on MSigDB Hallmark only (clusterProfiler GSEA / ORA / Reactome to follow) |
| `new_run_dir`, `write_run_record`, `render_run_report` (`R/report.R`, `inst/report/run_report.Rmd`) | `common.R` (`ensure_run_scaffold`: `00_input_freeze` pattern) | one folder per run with MD5 manifest, parameters, versions, tables, figures, HTML report |

Not ported (thesis-specific): stages 03c–03l (profile plots for the PRM panel, validation ranking,
pathway crosstalk, biology themes, evidence matrix, PRM selection).

## Demo data (`demo_data/`; the parquet and the FASTA are gitignored)

| file | source (read only) | how it was made |
|---|---|---|
| `report.parquet` (36 runs, 361 MB) | PhD `data_locked/inputs/diann_v231_output_msaccuracy_auto/report.parquet` (DIA-NN 2.3.1, 76 runs) | `scripts/make_demo_subset.R` (2026-11-04): Long-term (stage L) runs only; rows limited to the peptides MS-DAP keeps on the full 76-run import (best charge state per modified sequence, plain sequences mapped to one protein group); original row order preserved. |
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
