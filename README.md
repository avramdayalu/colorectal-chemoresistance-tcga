# Transcriptomic Signatures of Chemoresistance in Colorectal Cancer (TCGA-COAD/READ)

An independent bioinformatics project analyzing RNA-seq data from TCGA colon and rectal
adenocarcinoma patients (TCGA-COAD + TCGA-READ) to characterize gene expression differences
between chemoresistant and chemosensitive tumors, using differential expression analysis,
pathway enrichment, and — the central focus of this repository — a set of confound and
robustness checks that test whether the headline findings hold up under scrutiny.

This is a self-study project (R, DESeq2, GSEA), not a peer-reviewed study or clinical
research. It's presented here with its methods, results, and limitations stated as plainly
as possible, including analyses that complicated or partly undercut the initial findings.

## Why this repo looks the way it does

Most tutorial-style differential expression projects stop at "run DESeq2, run GSEA, report
significant pathways." This one doesn't, because an early result (a set of proliferation and
immune-related Hallmark pathways separating the two groups) turned out to have a real
confound sitting underneath it, once the underlying clinical treatment records were checked.
Chasing that confound down — quantifying it, testing whether the main result survives
removing it, and then testing whether *that* result was itself just a sample-size artifact —
is the actual substance of this project. The three-stage structure below (main analysis →
sensitivity analysis → permutation control) reflects that process in the order it happened.

## Data

- **Source:** TCGA-COAD (colon adenocarcinoma) and TCGA-READ (rectal adenocarcinoma),
  STAR-Counts RNA-seq quantification and clinical data, downloaded from the NCI Genomic
  Data Commons (GDC) via `TCGAbiolinks`.
- **Cohort size:** 87 primary tumor samples (19 chemoresistant / 68 chemosensitive).
- **Disease type restriction:** `cases.disease_type == "Adenomas and Adenocarcinomas"`
  (excludes rarer cystic/mucinous/serous histologic subtypes).
- Raw count matrices and clinical files are **not** committed to this repository (large,
  and redistribution of GDC data has its own terms — see [GDC Data Access
  policies](https://gdc.cancer.gov/access-data)). `download_tcga_colorectal_counts.R`
  documents exactly how the count matrix and clinical export were pulled, so the dataset
  is reproducible from GDC directly.

### How "chemoresistant" and "chemosensitive" are defined — read this before the results

This is the single most important methodological point in the project, and it wasn't
obvious from the outset, so it's stated here explicitly rather than left implicit.

Groups are defined by **treatment *outcome*, not by a specific drug or regimen**:

- **Chemoresistant (C):** `treatments.treatment_outcome == "Progressive Disease"`
- **Chemosensitive (S):** `treatments.treatment_outcome == "Complete Response"`

restricted to cases with chemotherapy treatment recorded
(`treatments.treatment_type %in% c("Chemotherapy", "Pharmaceutical Therapy, NOS")`) and a
single, unambiguous outcome label across all recorded lines of therapy. Cases with no
chemo-outcome record, conflicting outcomes across treatment lines, or an intermediate
outcome (Partial Response / Stable Disease / Treatment Ongoing) were excluded rather than
guessed at.

Querying the underlying GDC clinical records (`treatments.therapeutic_agents`) directly
shows the actual regimens are **heterogeneous, predominantly FOLFOX-based**
(fluorouracil + leucovorin + oxaliplatin, in various combinations), with a minority on
FOLFIRI, capecitabine monotherapy, or a targeted biologic added to the regimen
(bevacizumab, cetuximab, or ziv-aflibercept) — not one uniform drug or protocol. This
matters for interpreting the results: this is a study of a treatment-outcome-associated
expression signature across a heterogeneous real-world chemotherapy population, not a
controlled comparison of one specific drug's mechanism.

A second, equally important sampling caveat: TCGA primary tumor RNA-seq specimens are
mostly collected **before** treatment. In this cohort, `diagnoses.prior_treatment == "No"`
for 71/87 cases (82%); only 10/87 were sequenced after some prior treatment (plausibly
explained by neoadjuvant chemoradiation being standard of care for rectal but not colon
cancer). So the expression signatures reported here reflect a **baseline, pre-treatment
transcriptomic state associated with later treatment outcome** — not a real-time
pharmacodynamic response to drug exposure.

## Pipeline

All scripts are config-driven off `config/config.yaml` (sample ID lists, filtering
thresholds, DESeq2 contrast, GSEA/ORA parameters, and every input/output file path), run
in this order:

| Script | What it does |
|---|---|
| `src/data_loader.R` | Loads the raw count matrix, applies a low-count filter, builds the DESeq2 sample/condition table from the C/S case ID lists above. |
| `src/deseq2_pipeline.R` | Fits the DESeq2 negative-binomial model (design `~condition`, contrast S vs. C), PCA on variance-stabilized counts, extracts results, applies `apeglm` log2FC shrinkage, writes the full and filtered results tables. |
| `src/make_volcano_plots.R` | Volcano plots (all genes + labeled candidate genes). |
| `src/gsea_analysis.R` | Preranked GSEA (`fgsea`) against MSigDB Hallmark gene sets, genes ranked by the (unshrunk) Wald statistic. |
| `src/leading_edge_analysis.R` | Enrichment curves and per-patient leading-edge gene heatmaps for the top GSEA hits (a consistency check, not independent validation — see script header). |
| `src/ora_analysis.R` | Complementary Gene Ontology over-representation analysis (`clusterProfiler::enrichGO`, hypergeometric test, explicit background universe), run separately on genes up in S vs. up in C, as a cross-check against a different statistical method on the same data. |
| `src/sensitivity_analysis.R` | Re-runs DESeq2 + GSEA after excluding the subset of chemoresistant cases who received a targeted biologic on top of chemotherapy (see *Confound analysis* below). |
| `src/permutation_control.R` | Random-subsampling control: repeats the same exclusion procedure with 50 randomly chosen 8-case drops, to test whether the sensitivity-analysis result is specific to the biologic-treated cases or just a generic effect of a smaller sample. |

## Results

**Differential expression:** a small set of genes clear a combined significance and
effect-size threshold (`results/tables/candidate_genes_final.csv`; thresholds set in
`config.yaml`).

**GSEA (Hallmark, preranked):** 13 gene sets reach `padj ≤ 0.05`, splitting cleanly into
two clusters — proliferation/immune pathways (E2F_TARGETS, G2M_CHECKPOINT, MYC_TARGETS_V1,
INTERFERON_GAMMA_RESPONSE, INTERFERON_ALPHA_RESPONSE, INFLAMMATORY_RESPONSE,
ALLOGRAFT_REJECTION, COMPLEMENT) enriched in the chemoresistant group, and a
metabolic/lipid cluster (CHOLESTEROL_HOMEOSTASIS and related pathways) enriched in the
chemosensitive group. Full table: `results/tables/gsea_hallmark_results.csv`.

**ORA cross-check:** genes up in the chemosensitive group are independently enriched for
lipid transport/localization Gene Ontology terms — consistent with the GSEA
CHOLESTEROL_HOMEOSTASIS result, via a different statistical test on the same expression
data.

**PCA:** whole-transcriptome PCA on variance-stabilized counts does **not** show a clean
separation between chemoresistant and chemosensitive samples — the pathway-level signal
above is real in the statistical sense (padj-adjusted), but it is not a dominant, globally
obvious axis of variation in this dataset. Worth stating plainly rather than glossing over.

## Confound analysis: how much of this is really about treatment outcome, and how much is a treatment-exposure artifact?

The heterogeneity described above raises an obvious question: does the C-vs-S signature
partly reflect *which drugs were used* rather than the biology of resistance itself? One
concrete, testable version of this: targeted biologics (bevacizumab, cetuximab,
ziv-aflibercept) are typically added only after a tumor progresses on standard
chemotherapy — so biologic use should be systematically associated with the chemoresistant
group, and could plausibly be driving part of the immune/inflammatory pathway signal
(these drugs affect the tumor microenvironment) rather than reflecting baseline tumor
biology.

Querying `treatments.therapeutic_agents` directly confirms the association: **8/19 (42.1%)
chemoresistant cases** received a biologic, versus **4/68 (5.9%) chemosensitive cases**
(Fisher's exact test, OR = 11.64, p = 0.0004).

### Sensitivity analysis

`src/sensitivity_analysis.R` re-runs the full DESeq2 + GSEA pipeline after excluding the 8
biologic-treated chemoresistant cases (n=79: 11 C vs. 68 S). The result was not a simple
confirmation or refutation — it was a real shift in both directions:

| Pathway | padj, n=87 (original) | padj, n=79 (biologics excluded) |
|---|---|---|
| E2F_TARGETS | 4.8 × 10⁻¹⁷ | 0.46 (lost significance) |
| MYC_TARGETS_V1 | 0.0032 | 1.00 (lost significance, NES sign flipped) |
| G2M_CHECKPOINT | 1.8 × 10⁻¹⁵ | 0.0036 (weakened, still significant) |
| INTERFERON_GAMMA_RESPONSE | 1.2 × 10⁻⁷ | 3.2 × 10⁻²⁴ (strengthened) |
| INFLAMMATORY_RESPONSE | 0.048 | 2.9 × 10⁻¹⁴ (strengthened) |
| ALLOGRAFT_REJECTION | 0.005 | 1.3 × 10⁻¹⁷ (strengthened) |
| COMPLEMENT | 0.02 | 8.9 × 10⁻¹⁰ (strengthened) |
| CHOLESTEROL_HOMEOSTASIS | 0.005 | 0.052 (weakened to borderline) |

Full table: `results/tables/gsea_sensitivity_comparison.csv`. Two proliferation pathways
collapsed; four immune pathways got *more* significant, not less. Neither the simple
"proliferation is real baseline biology, immune is a treatment artifact" story nor its
reverse fits this cleanly — and going from 19 to 11 chemoresistant cases is itself a large
enough change in statistical power that some of this shift could just be noise. That
ambiguity is exactly what the next step tests directly, rather than asserting an
interpretation from a single before/after comparison.

### Permutation control

`src/permutation_control.R` repeats the exact same exclusion procedure 50 times, each time
dropping a different **random** set of 8 chemoresistant cases (not the biologic-treated
ones specifically), and compares the real biologics-exclusion result against that empirical
null distribution. This is the actual answer to "is this specific to the biologic-treated
cases, or just what happens when you shrink the group?":

- **MYC_TARGETS_V1's collapse is generic, not biologics-specific.** Exactly 50% of random
  8-case drops also lose MYC significance, with the sign flipping constantly — the real
  result is unremarkable against this backdrop. This rules out "the biologic-treated cases
  specifically carried the MYC signal."
- **E2F_TARGETS is a weaker, suggestive case.** Only 12% of random drops (6/50) reproduce
  a similar loss of significance — rarer than average, but not below a conventional
  significance threshold on its own.
- **G2M_CHECKPOINT is robust regardless of which cases are dropped** (94% of random draws
  stay significant, matching the real result) — never seriously in question.
- **Three of the five immune pathways strengthen generically** (interferon-alpha response,
  allograft rejection, complement): the real result falls comfortably inside the range
  produced by random draws, so their strengthening looks like a function of shrinking and
  homogenizing the smaller group, not something specific to biologics.
- **Two pathways are genuine outliers: INTERFERON_GAMMA_RESPONSE and
  INFLAMMATORY_RESPONSE.** The real biologics-exclusion NES for both is more extreme than
  *every one* of the 50 random draws (empirical one-sided p < 1/50 ≈ 0.02 for each). This
  is real, specific evidence that something about the biologic-treated cases — not just
  "having fewer cases" — is contributing to these two pathways' signal.

Full per-permutation results: `results/permutation_control/tables/permutation_gsea_results.csv`;
summary: `results/permutation_control/tables/permutation_summary.csv`.

**Bottom line:** the honest read of this whole confound-analysis thread is that most of
what looked like a clean "immune vs. proliferation" resistance signature does not survive
scrutiny as originally framed. The proliferation signal is partly a treatment-exposure
artifact (MYC) and partly underpowered rather than definitively real or fake (E2F). Most of
the immune-pathway strengthening after removing biologic-treated cases is a generic
small-sample effect, not evidence of an "unmasked" homogeneous immune subtype — except for
interferon-gamma response and inflammatory response specifically, where there's a real,
quantifiable, biologics-associated signal worth further investigation.

## Limitations

- **Single retrospective cohort, no independent validation set.** All confound checks
  above are internal to this same 87-sample cohort.
- **Associative, not causal or mechanistic.** Nothing here establishes *why* the
  interferon-gamma/inflammatory association with biologic-treated resistant cases exists —
  only that it's statistically real and specific in this cohort.
- **MSI (microsatellite instability) status was not available** in the downloaded GDC
  clinical export and was not checked as a possible confound, despite MSI-high colorectal
  tumors being independently associated with strong baseline interferon-gamma/immune
  signatures. A concrete, well-scoped next step if this project is extended.
- **Tumor purity / stromal-immune infiltration (e.g. ESTIMATE scores)** was not computed;
  bulk RNA-seq immune pathway signals can partly reflect microenvironment composition
  rather than tumor-intrinsic biology.
- **n=11 in the sensitivity/permutation analyses is small**, even though the permutation
  control's use of an empirical null (rather than trusting a single padj value) is the
  right way to handle that, not a way around it.

## Reproducing this analysis

Requires R with `DESeq2`, `fgsea`, `msigdbr`, `clusterProfiler`, `org.Hs.eg.db`,
`enrichplot`, `pheatmap`, `here`, `yaml`, `ggplot2` (Bioconductor + CRAN).

```r
source("src/deseq2_pipeline.R"); pipeline_out <- run_pipeline()
source("src/gsea_analysis.R");   gsea_out <- run_gsea(pipeline_out)
source("src/leading_edge_analysis.R"); run_leading_edge_analysis(gsea_out)
source("src/ora_analysis.R");    run_ora(pipeline_out)
source("src/sensitivity_analysis.R")   # after the above has run at least once
source("src/permutation_control.R")    # takes ~35 min for n_perm=50; edit that value to change
```

Raw count/clinical data are not included (see *Data* above) — `download_tcga_colorectal_counts.R`
documents how they were retrieved from GDC.

## Methods and citations

- Love MI, Huber W, Anders S (2014). Moderated estimation of fold change and dispersion
  for RNA-seq data with DESeq2. *Genome Biology*, 15(12), 550.
- Zhu A, Ibrahim JG, Love MI (2019). Heavy-tailed prior distributions for sequence count
  data: removing the noise and preserving large differences. *Bioinformatics*, 35(12),
  2084–2092. (apeglm shrinkage)
- Subramanian A, Tamayo P, Mootha VK, et al. (2005). Gene set enrichment analysis: a
  knowledge-based approach for interpreting genome-wide expression profiles. *PNAS*,
  102(43), 15545–15550.
- Korotkevich G, Sukhov V, Budin N, et al. (2021). Fast gene set enrichment analysis.
  *bioRxiv* 060012 (preprint, fgsea implementation).
- Liberzon A, Birger C, Thorvaldsdottir H, et al. (2015). The Molecular Signatures
  Database Hallmark gene set collection. *Cell Systems*, 1(6), 417–425.
- Yu G, Wang LG, Han Y, He QY (2012). clusterProfiler: an R package for comparing
  biological themes among gene clusters. *OMICS*, 16(5), 284–287; updated in Wu T, Hu E,
  Xu S, et al. (2021). *The Innovation*, 2(3), 100141.
- Ashburner M, Ball CA, Blake JA, et al. (2000). Gene Ontology: tool for the unification
  of biology. *Nature Genetics*, 25(1), 25–29.
- Kolde R. pheatmap (CRAN package, cited as software; no dedicated methods paper).
- TCGA data via the NCI Genomic Data Commons: https://gdc.cancer.gov

## Repository structure

```
config/config.yaml          # all paths, thresholds, design/contrast, GSEA/ORA parameters
src/                         # analysis scripts, run in the order listed above
data/raw/                    # (not committed) raw count matrix + gene index, from GDC
results/
  tables/                    # DESeq2 results, GSEA/ORA tables, sensitivity comparison
  figures/                   # PCA, volcano plots, GSEA plot, enrichment curves, heatmaps
  sensitivity_no_biologics/  # results from the biologics-excluded (n=79) rerun
  permutation_control/       # results from the 50-permutation random-drop control
```
