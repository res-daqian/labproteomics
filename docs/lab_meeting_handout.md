# Lab Proteomics Workbench: what it is, what it does, what it needs from you

*Working name. Handout for the lab meeting of 2026-11-06. Written in plain English; questions to Daqian.*

## 1. Why we built it

Every proteomics project in the lab goes through the same steps: import the search-engine
output, check the runs, normalise, compare groups, interpret, and keep a record of what was
done. Today these steps live in scripts that one person can run. The workbench packages the
same analysis chain as a local application, so that anyone in the lab can run it on their own
data without writing code, with the same defaults, and with every analysis written to a folder
that can be opened months later.

Three ideas shape it:

- **Quality control before statistics.** Runs are judged on three layers (design, single run,
  whole cohort) and exclusions are decided before any comparison is made.
- **Every choice is visible.** Each threshold is a parameter in the interface and in the record.
  Nothing is hidden in code.
- **Every hit is described, not only counted.** A hit comes with its exact rule, its per-run
  values and a stability check.

## 2. What goes in, what comes out

| In | Out |
|---|---|
| The main report of the search engine (DIA-NN `report.parquet` or `report.tsv`; Spectronaut and FragPipe outputs are connected but not yet tested with data) | A run record: one folder per analysis with the input checksum, the annotation, all parameters, package versions, tables, figures (PDF and PNG), an HTML report and a PDF QC report |
| Optional: the search engine's per-run statistics (`report.stats.tsv`) and log (`report.log.txt`) next to the report | Upload tables for IPA, iPathwayGuide, STRING and Cytoscape |
| A sample annotation table (CSV or Excel): run name, group, subject, batch, run order | An mzQC file with the per-run metrics (HUPO-PSI standard) |

The application runs on your own computer, inside R. Nothing is uploaded anywhere.

## 3. The six pages

| Page | Question it answers | What you decide |
|---|---|---|
| 1 Data | Is the design analysable? Which runs belong to which group, subject and batch? | Attach the annotation; read the design check (batch against group, pairing, run order, pooled QC) |
| 2 Quality control | Which runs are weak, contaminated or unusual? | Exclude runs, from evidence, before any comparison |
| 3 Differential analysis | Which proteins differ between two groups? | The contrasts, paired or unpaired; the gate (fold change with p, q or limma treat) |
| 4 Rescue (only with three or more groups) | Which disease-perturbed proteins does a treatment bring back? | The disease contrast and the treatment contrast |
| 5 Enrichment | Which gene sets move? | GSEA on all proteins, or ORA of the hits against the quantified background |
| 6 Run record | What exactly was done? | Label and notes; everything else is written automatically |

## 4. Quality control in three layers

- **Design (page 1).** Before any number is computed: are batches crossed with the groups or
  confounded with them; can the contrasts be paired; is the acquisition order mixed; are pooled
  QC injections present. A confounded contrast is still analysed, but it is labelled supportive.
- **Single run (page 2, scorecard).** The search engine's own per-run statistics: identifications,
  signal, peak width, mass error, retention-time deviation, normalisation instability, missed
  cleavages, peptide length and charge. Each run is compared with the other runs of its batch
  (median and MAD); batch-wide differences are reported separately as systematic effects.
- **Composition (page 2, marker panel).** The share of the signal carried by marker proteins of
  the tissue itself, of neighbouring tissues (for example lens or RPE in a retina sample), of
  blood and plasma, and of handling contaminants (skin keratins). A Y-chromosome set checks the
  recorded sex of each sample.
- **Cohort (page 2, after normalisation).** Missing values, intensity centring, PCA, sample
  correlation, acquisition drift, batch structure and within-subject consistency, with MAD-based
  flags. A run with two or more flags is a candidate for exclusion.

The rule is simple: look at all of this first, write down which runs you exclude and why, and
only then run the comparisons. The record keeps the decision.

## 5. How a hit is defined

The default gate is the one used in the retina atropine study:

```
|log2 fold change| >= log2(1.2) (= 0.263) AND p < 0.05
```

with DEqMS as the model (peptide-count-aware moderated t-test on MaxLFQ protein intensities after
vwmb normalisation). The gate can be switched to the Benjamini-Hochberg q-value, or to limma
treat, which tests the fold-change threshold directly. Three extra descriptions come with every
contrast and never change the hit list:

- **Fragile.** The hit disappears when one subject (or one run) is left out and the model is
  refitted. With small groups this is common; it marks the evidence level, not an error.
- **Detection-only.** Proteins seen in one group and almost never in the other have no reliable
  fold change. They are listed separately with a Fisher exact test.
- **limma treat.** How many hits survive a test of the fold-change threshold itself.

## 6. Rescue and enrichment

- **Rescue** needs a disease contrast (disease versus control) and a treatment contrast (treated
  disease versus disease). Tier 1: the treatment reverses the disease change with its own
  significance. Tier 2: the treatment is not significant but the protein is back within the
  fold-change band. The definitions are printed next to the counts. A drug-only contrast adds a
  flag and never filters.
- **Enrichment** offers GSEA (all proteins ranked by the model t-statistic) and ORA (hits against
  the proteins quantified in that contrast, never against the whole genome). Collections: MSigDB
  Hallmark, GO, Reactome, KEGG MEDICUS and cell-type signatures, offline. Non-human data use an
  ortholog table to human symbols.

## 7. What it does not do yet

- Spectronaut and FragPipe inputs are connected to the MS-DAP importers but have not been tested
  with real files; the run scorecard reads DIA-NN statistics only.
- The leave-one-out check refits the statistical model; it does not repeat normalisation.
- No phenotype association, no WGCNA, no second statistical engine yet.
- Computation is synchronous: a progress bar, about half a minute for QC and for each contrast on
  a few dozen runs.

## 8. What we need from you

- **Your annotation table.** For every run: a short name, the group, the subject (animal or
  donor; needed for paired designs), the batch (digestion date or processing batch) and the
  acquisition order. The Data page gives you a template with the run names filled in.
- **Pooled QC injections in future experiments.** A pooled sample injected several times across
  the sequence is the only direct measure of technical precision and drift. Old data sets without
  it still work; the workbench falls back to within-subject consistency and the drift plot.
- **Marker lists for your tissue.** The built-in sets cover retina, RPE and choroid, sclera, tears
  and cultured endothelial cells as a start. Tell us which proteins mark your tissue and its usual
  contaminants; a custom list is a two-column CSV.
- **Your data types and comparisons.** DIA-NN output is tested; Spectronaut users are the next
  test case.

## 9. How to get it

Install R and RStudio, open the workbench folder, run `source("scripts/install_dependencies.R")`
once, then install the package (README, section "Install"). The demo data set is public (PRIDE
PXD079891, human endothelial cells, 16 runs) and is downloaded with one command. Start the
application by double-clicking `launch_mac.command` or `launch_windows.bat`.

## 10. Questions people ask

- **Is my data uploaded?** No. Everything runs inside R on your computer. The browser is only the
  interface.
- **Why p instead of q by default?** With paired designs of 6 to 10 subjects the q-value is very
  strict. Both columns are always shown, and the gate can be switched in one click.
- **How is pairing handled?** The subject becomes a blocking factor in the model; only subjects
  with a run in both groups are used.
- **Can I change the defaults?** Yes, every one of them, and the record shows what you used.
- **What if my data set has only two groups?** Everything works; the Rescue page simply does not
  appear.
