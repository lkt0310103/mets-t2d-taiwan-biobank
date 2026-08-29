#!/bin/bash
# =============================================================================
# Single-trait genome-wide association analyses
#
# Association testing is performed separately in the training (60%) and
# testing (40%) sets for each obesity- and metabolic syndrome-related
# phenotype, using an additive model in PLINK 2.0.
#
# All continuous traits were standardised (Z score transformed) before
# analysis. Covariates: sex, age, smoking status, alcohol use, betel nut
# chewing, physical activity, genotyping batch and genetic principal
# components.
#
# Input : qc{6,4}-final.{bed,bim,fam}   -- post-QC genotypes (see 01_qc/)
#         <trait><set>.tab             -- FID IID phenotype (column 3)
#         covar<set>.txt               -- covariate file
# Output: result_<trait><set>.<pheno>.glm.linear (or .glm.logistic)
#
# Software: PLINK v2.0.0-a.6.8  
# =============================================================================

set -euo pipefail

# ---------------------------------------------------------------------------
# 0. Configuration
# ---------------------------------------------------------------------------
WORKDIR="./data"
THREADS=48
MEMORY=100000          # MB
PLINK2="plink2"

# Covariates. The number of principal components retained was decided from a
# the eigenvalues thrshold produced in 01_qc/ (see manuscript ESM Methods).
#

COVARS_TRAIN="age,genotype_BATCH,PC1-PC6,DRINKING,PHY,SMOKINGG,S,nut"
COVARS_TEST="age,genotype_BATCH,PC1-PC4,DRINKING,PHY,SMOKINGG,S,nut"

# Continuous phenotypes

TRAITS_CONT="bmi bf wc whr sbp dbp hdl ldl tg fg hba1c"

cd "${WORKDIR}"

# ---------------------------------------------------------------------------
# 1. Continuous traits
# ---------------------------------------------------------------------------
for SET in 6 4; do

  if [ "${SET}" = "6" ]; then
    COVARS="${COVARS_TRAIN}"
  else
    COVARS="${COVARS_TEST}"
  fi

  for TRAIT in ${TRAITS_CONT}; do

    # Restrict to individuals with a non-missing value for this trait
    ${PLINK2} \
      --bfile "qc${SET}-final" \
      --threads ${THREADS} \
      --keep "${TRAIT}${SET}.tab" \
      --make-bed \
      --out "${TRAIT}${SET}"

    # Attach the phenotype (column 3 of the .tab file)
    ${PLINK2} \
      --bfile "${TRAIT}${SET}" \
      --threads ${THREADS} \
      --pheno "${TRAIT}${SET}.tab" \
      --pheno-col-nums 3 \
      --make-bed \
      --out "${TRAIT}${SET}_gwas"

    # Additive linear model
    ${PLINK2} \
      --bfile "${TRAIT}${SET}_gwas" \
      --threads ${THREADS} \
      --memory ${MEMORY} \
      --glm hide-covar \
      --covar "covar${SET}.txt" \
      --covar-name "${COVARS}" \
      --covar-variance-standardize \
      --out "result_${TRAIT}${SET}"

  done
done

# ---------------------------------------------------------------------------
# 2. Type 2 diabetes (binary trait, logistic model)
#
# ---------------------------------------------------------------------------
# for SET in 6 4; do
#   ${PLINK2} \
#     --bfile "qc${SET}-final" \
#     --threads ${THREADS} \
#     --memory ${MEMORY} \
#     --pheno "t2d${SET}.tab" --pheno-col-nums 3 \
#     --glm firth-fallback hide-covar \
#     --covar "covar${SET}.txt" \
#     --covar-name "${COVARS}" \
#     --covar-variance-standardize \
#     --out "result_t2d${SET}"
# done

echo "GWAS complete."
