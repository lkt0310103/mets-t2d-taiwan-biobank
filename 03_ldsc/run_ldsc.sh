#!/bin/bash
# =============================================================================
# Linkage disequilibrium score regression (LDSC)
#
# Estimates SNP heritability, the LDSC intercept and pairwise genetic
# correlations from the training-set GWAS summary statistics.
#
#   - Univariate LDSC : all 12 phenotypes (11 continuous traits + type 2 diabetes)
#   - Bivariate LDSC  : the 55 pairs among the 11 continuous traits
#   - Type 2 diabetes heritability is additionally converted to the liability
#     scale (6262 cases, 57,512 control participants; population prevalence
#     10.68%)
#
# Software: LDSC v1.0.1 (https://github.com/bulik/ldsc), Python 2.7
# Reference: 1000 Genomes Project Phase 3 East Asian (EAS) LD scores
#
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# 0. Configuration
# ---------------------------------------------------------------------------
LDSC_DIR="./ldsc"                       # cloned LDSC repository
SUMSTATS_DIR="./data/sumstats"          # PLINK 2.0 --glm output, training set
MUNGED_DIR="./results/munged"
OUT_DIR="./results/ldsc"

EAS_LD="./reference/eas_ldscores/"      # pre-computed 1000G EAS LD scores
HAPMAP_SNPS="./reference/w_hm3.snplist" # SNP list used for merge-alleles

mkdir -p "${MUNGED_DIR}" "${OUT_DIR}"

# 11 continuous traits + type 2 diabetes
TRAITS="bmi bf wc whr sbp dbp hdl ldl tg fg hba1c t2d"

# ---------------------------------------------------------------------------
# 1. Pre-process summary statistics
#
# munge_sumstats.py removes variants without allele frequency information,
# non-biallelic SNPs, indels, MHC-region variants and SNPs with MAF < 1%,
# and restricts to the reference SNP list. Variants with an imputation
# quality (INFO) score < 0.8 were excluded upstream.
# ---------------------------------------------------------------------------
for TRAIT in ${TRAITS}; do


  python "${LDSC_DIR}/munge_sumstats.py" \
    --sumstats "${SUMSTATS_DIR}/result_${TRAIT}6.txt" \
    --snp      ID \
    --a1       A1 \
    --a2       OMITTED \
    --frq      A1_FREQ \
    --p        P \
    --signed-sumstats BETA,0 \
    --N-col    OBS_CT \
    --maf-min  0.01 \
    --info-min 0.8 \
    --merge-alleles "${HAPMAP_SNPS}" \
    --chunksize 500000 \
    --out "${MUNGED_DIR}/${TRAIT}"

done

# ---------------------------------------------------------------------------
# 2. Univariate LDSC -- SNP heritability and intercept
# ---------------------------------------------------------------------------
for TRAIT in ${TRAITS}; do

  python "${LDSC_DIR}/ldsc.py" \
    --h2 "${MUNGED_DIR}/${TRAIT}.sumstats.gz" \
    --ref-ld-chr "${EAS_LD}" \
    --w-ld-chr   "${EAS_LD}" \
    --out "${OUT_DIR}/h2_${TRAIT}"

done

# Type 2 diabetes: convert observed-scale heritability to the liability scale
python "${LDSC_DIR}/ldsc.py" \
  --h2 "${MUNGED_DIR}/t2d.sumstats.gz" \
  --ref-ld-chr "${EAS_LD}" \
  --w-ld-chr   "${EAS_LD}" \
  --samp-prev 0.0982 \
  --pop-prev  0.1068 \
  --out "${OUT_DIR}/h2_t2d_liability"
# samp-prev = 6262 / (6262 + 57512) = 0.0982
# pop-prev  = 10.68% (Taiwan diabetes atlas 2024)

# ---------------------------------------------------------------------------
# 3. Bivariate LDSC -- genetic correlations
#
# All 55 pairs among the 11 continuous traits. Type 2 diabetes is excluded
# from the genetic correlation and downstream pleiotropy analyses.
#
# Significance threshold: Bonferroni-corrected p < 0.05/55 = 9.09e-4
# ---------------------------------------------------------------------------
CONT_TRAITS=(bmi bf wc whr sbp dbp hdl ldl tg fg hba1c)

for i in "${!CONT_TRAITS[@]}"; do
  for j in "${!CONT_TRAITS[@]}"; do
    if [ "$i" -lt "$j" ]; then
      T1="${CONT_TRAITS[$i]}"
      T2="${CONT_TRAITS[$j]}"

      python "${LDSC_DIR}/ldsc.py" \
        --rg "${MUNGED_DIR}/${T1}.sumstats.gz,${MUNGED_DIR}/${T2}.sumstats.gz" \
        --ref-ld-chr "${EAS_LD}" \
        --w-ld-chr   "${EAS_LD}" \
        --out "${OUT_DIR}/rg_${T1}_${T2}"

    fi
  done
done

echo "LDSC complete. Results in ${OUT_DIR}"
