# Genetic architecture of metabolic syndrome and type 2 diabetes risk prediction in the Taiwan Biobank

Analysis code accompanying:

> Li K-T, Tai T-H, Wang T-N. *Dissecting the genetic architecture of metabolic
> syndrome to improve type 2 diabetes risk prediction: factor-based GWAS and
> Mendelian randomisation.* 

---

## Data availability

Individual-level genotype and phenotype data were obtained from the
**Taiwan Biobank** and are **not redistributable**. Access can be applied for
through the Taiwan Biobank (https://www.twbiobank.org.tw/) subject to approval
by its Ethics and Governance Council. No participant-level data are contained
in this repository.

GWAS summary statistics generated in this study are available on request from
the corresponding author.

Ethics approval: Institutional Review Board of Kaohsiung Medical University
Hospital, KMUHIRB-G(I)-20180009.

---

## Analysis pipeline

Scripts are numbered in the order they were run. Throughout, the suffix `6`
denotes the training set (60% of the cohort) and `4` the testing set (40%).

| Directory | Contents |
|---|---|
| `01_qc/` | Sample- and variant-level quality control; principal component analysis |
| `02_gwas/` | Single-trait genome-wide association analyses |
| `03_ldsc/` | SNP heritability, LDSC intercept and pairwise genetic correlations |
| `04_factor_analysis/` | Exploratory factor analysis (genetic correlation matrix) and second-order confirmatory factor analysis (testing set) |
| `05_factor_gwas/` | Factor-based GWAS of the latent MetS factor using Genomic SEM |
| `06_pleiotropy/` | SNP- and gene-level PLACO+ across 55 trait pairs; MAGMA gene-based analysis |
| `07_mendelian_randomisation/` | Bidirectional two-sample Mendelian randomisation |
| `08_prs/` | Polygenic risk score construction and validation |

### Order of execution

```
01_qc/run_qc.sh
02_gwas/run_gwas.sh
03_ldsc/run_ldsc.sh
04_factor_analysis/efa.R
04_factor_analysis/cfa.R
05_factor_gwas/genomic_sem_gwas.R
06_pleiotropy/placo_snp_level.R
06_pleiotropy/magma_gene_level.sh
06_pleiotropy/placo_gene_level.R
07_mendelian_randomisation/bidirectional_mr.R
08_prs/prs_pipeline.sh
08_prs/prs_evaluation.R
```

---

## Software


| Software | Version | Notes |
|---|---|---|
| PLINK | v1.9.0-b.7.7 | Sample QC, LD pruning, relatedness, clumping |
| PLINK | v2.0.0-a.6.8 | Dataset handling, association testing, scoring |
| LDSC | 1.0.1 | https://github.com/bulik/ldsc |
| R | 4.5.1 | |
| `psych` | 2.6.5 | Exploratory factor analysis |
| `lavaan` | 0.7-2 | Confirmatory factor analysis |
| `GenomicSEM` | 0.0.5c | Factor-based GWAS |
| PLACO+ | v0.2.0. | Park & Ray 2026 |
| MAGMA | v1.10 | Gene-based analysis |
| `TwoSampleMR` | 0.7.9 | Mendelian randomisation |
| `pROC` | 1.19.1 | AUC, DeLong's test |
| `fmsb` | 0.7.6 | Nagelkerke pseudo-R² |
| FUMA | 2.1.5. | Web application; see below |
| SAS | 9.4 | Baseline characteristics |

### Reference data

- 1000 Genomes Project Phase 3 East Asian (EAS) reference panel — used as the
  linkage disequilibrium reference throughout
- Pre-computed EAS LD scores — LDSC and Genomic SEM
- HapMap3 SNP list — Genomic SEM harmonisation
- NCBI Build 38 gene locations, 18,469 protein-coding genes — MAGMA
- Ensembl v.110 — FUMA background gene set
- GTEx v8, 54 specific tissue types — FUMA tissue expression profiling

---

## FUMA settings

Functional annotation and gene-set enrichment were performed with the
GENE2FUNC module of FUMA (https://fuma.ctglab.nl/), which is a web
application and therefore has no accompanying script. Settings used:

- Background gene set: Ensembl v.110 (protein-coding genes, long non-coding
  RNAs, non-coding RNAs, processed transcripts, pseudogenes, and
  immunoglobulin and T-cell receptor genes)
- Gene expression: GTEx v8, 54 specific tissue types, log₂-transformed
  average expression
- Differentially expressed gene enrichment: two-sided tests,
  Bonferroni-corrected p < 0.05
- Gene-set enrichment: Gene Ontology (BP, CC, MF), KEGG and Reactome;
  minimum overlap of two genes; Benjamini–Hochberg FDR-adjusted p < 0.05

Two gene sets were submitted separately: the 13 genes reaching gene-level
significance in the MAGMA analysis of the factor GWAS, and the 626 pleiotropic
genes identified at FDR < 0.05 in the gene-level PLACO+ analysis.

---

## Notes on reproducibility

- File paths are set in a configuration block at the top of each script; no
  absolute paths are hard-coded.
- Scripts marked `VERIFY` contain steps that should be checked against the
  input file formats in your own environment.
- Continuous net reclassification improvement was estimated with 50 bootstrap
  iterations under a fixed random seed.

---

## Licence

MIT — see `LICENSE`.
