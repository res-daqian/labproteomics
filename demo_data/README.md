# Demo data: PRIDE PXD079891 (public, CC0)

The workbench demo uses a public DIA-NN data set, not data from this laboratory.

| Item | Value |
|---|---|
| PRIDE accession | [PXD079891](https://www.ebi.ac.uk/pride/archive/projects/PXD079891) |
| Publication | Zhang SY et al., *The Journal of Physiology* 2026, DOI [10.1113/JP292141](https://doi.org/10.1113/JP292141) |
| Samples | Human umbilical vein endothelial cells (HUVEC) from 8 donors, 4 female (SZF1 to SZF4) and 4 male (SZM1 to SZM4) |
| Design | Each donor preparation transfected with a scrambled antisense oligonucleotide (group `Control`) or an ASO against the lncRNA MANCR (group `MANCR_KD`); 16 runs, paired by donor (`subject_id`) |
| Instrument and search | Orbitrap Exploris 480, DIA; DIA-NN 1.9, library-free search against human Swiss-Prot (2024-11-11), 1 percent FDR |
| Licence | PRIDE data are released under CC0 |

## Files

| File | Source | Size (bytes) |
|---|---|---|
| `report.parquet` | PRIDE, main DIA-NN report (gitignored; `scripts/get_demo_data.R` downloads it) | 164024858 |
| `report.stats.tsv` | PRIDE, DIA-NN per-run statistics (used by the run scorecard) | 3645 |
| `report.log.txt` | PRIDE, DIA-NN log (version and settings shown on the Data page) | 44205 |
| `PXD079891_community_annotated.sdrf.tsv` | PRIDE, sample metadata in SDRF format | 13838 |
| `demo_sample_annotation.csv` | built from the SDRF by `scripts/get_demo_data.R` | small |
| `demo_settings.json` | tissue and species used by the demo button | small |
| `UP000005640_9606.fasta` | UniProt human reference proteome, optional (gitignored); gives protein descriptions | about 13 MB |

## Annotation columns

`sample_id` is the DIA-NN run name (`MtoZ-MZJL0304202603-01` to `-16`); `sample_name` is the
authors' label (for example `SZF1NC` = donor SZF1, control). `subject_id` is the donor, so the
knockdown contrast is paired. `batch` holds the donor sex: the female and male preparations were
acquired as two consecutive blocks (runs 1 to 8 and 9 to 16), so sex is the only grouping of runs
that can act like a batch here; the same value is kept in `sex` for the sex check of the marker
panel. No run is excluded in the shipped annotation. The authors excluded donor SZF3 (runs 05 and
06) from their differential analysis as a PCA outlier; the QC page lets you take that decision
yourself and records it.

## Rebuild

```
Rscript scripts/get_demo_data.R            # download from PRIDE
Rscript scripts/get_demo_data.R <folder>   # or copy from a folder that already holds the files
```
