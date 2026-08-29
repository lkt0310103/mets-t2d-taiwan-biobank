#!/bin/bash
# =============================================================================
# Polygenic risk score construction
#
# Three scores are built with a clumping-and-thresholding (C+T) approach.
# Clumping is performed in the training set (r2 < 0.1, 250 kb window); scores
# are computed in the testing set.
#
#   1. MetS-weighted factor PRS : SNPs from the latent MetS factor GWAS,
#                                 weighted by factor GWAS effect sizes
#   2. T2D-weighted factor PRS  : the same SNP set, reweighted by training-set
#                                 type 2 diabetes GWAS effect sizes
#   3. PLACO-based PRS          : the 31,966 SNPs meeting FDR < 0.05 in the
#                                 SNP-level PLACO+ analysis, clumped to 1147
#                                 independent variants, weighted by
#                                 training-set type 2 diabetes effect sizes
#
# Software: PLINK v1.90b6.x (clumping), PLINK v2.00a3.x (scoring)
# =============================================================================

set -euo pipefail

WORKDIR="./data"
THREADS=48
PLINK1="plink"
PLINK2="plink2"

TRAIN="qc6-final"
TEST="qc4-final"

FACTOR_GWAS="results/factor_gwas/factor_gwas_all.txt"     # SNP, P, BETA
T2D_GWAS="data/sumstats/result_t2d6.txt"
PLACO_SNPS="results/pleiotropy/snp_level/pleiotropic_snps_fdr05.txt"

OUT_DIR="results/prs"
mkdir -p "${OUT_DIR}"

THRESHOLDS="5e-8 1e-5 1e-3 0.01 0.05 0.1"

cd "${WORKDIR}"

# ---------------------------------------------------------------------------
# 1. Factor-derived scores: clump at each p value threshold
# ---------------------------------------------------------------------------
for P in ${THRESHOLDS}; do

  ${PLINK1} \
    --bfile "${TRAIN}" \
    --threads ${THREADS} \
    --clump "${FACTOR_GWAS}" \
    --clump-snp-field SNP \
    --clump-field P \
    --clump-p1 "${P}" \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --out "${OUT_DIR}/factor_clump_${P}"

  awk 'NR>1 && $3 != "" {print $3}' "${OUT_DIR}/factor_clump_${P}.clumped" \
    > "${OUT_DIR}/factor_snps_${P}.txt"

  echo "Threshold ${P}: $(wc -l < "${OUT_DIR}/factor_snps_${P}.txt") SNPs"
done

# ---------------------------------------------------------------------------
# 2. PLACO-based score: clump the full FDR < 0.05 SNP set
#
# The whole FDR-based set is clumped with the same LD parameters, yielding
# 1147 independent variants.
# ---------------------------------------------------------------------------
# Build a p value file restricted to the PLACO SNPs, using the type 2 diabetes
# GWAS p values to rank variants during clumping.
# VERIFY: which p value column was used to rank SNPs when clumping the PLACO
#         set (type 2 diabetes GWAS p, or the PLACO p value).
awk 'NR==FNR {keep[$1]; next} FNR==1 || ($3 in keep)' \
    "${PLACO_SNPS}" "${T2D_GWAS}" > "${OUT_DIR}/placo_for_clump.txt"

${PLINK1} \
  --bfile "${TRAIN}" \
  --threads ${THREADS} \
  --clump "${OUT_DIR}/placo_for_clump.txt" \
  --clump-snp-field ID \
  --clump-field P \
  --clump-p1 1 \
  --clump-r2 0.1 \
  --clump-kb 250 \
  --out "${OUT_DIR}/placo_clump"

awk 'NR>1 && $3 != "" {print $3}' "${OUT_DIR}/placo_clump.clumped" \
  > "${OUT_DIR}/placo_snps.txt"

echo "PLACO score: $(wc -l < "${OUT_DIR}/placo_snps.txt") independent SNPs"

# ---------------------------------------------------------------------------
# 3. Score in the testing set
#
# --score expects: SNP ID, effect allele, effect size
# ---------------------------------------------------------------------------
for P in ${THRESHOLDS}; do

  # MetS weights (factor GWAS effect sizes)
  ${PLINK2} \
    --bfile "${TEST}" \
    --threads ${THREADS} \
    --extract "${OUT_DIR}/factor_snps_${P}.txt" \
    --score "${FACTOR_GWAS}" 1 2 3 header cols=+scoresums \
    --out "${OUT_DIR}/prs_factor_metsweight_${P}"

  # T2D weights (type 2 diabetes GWAS effect sizes, same SNP set)
  ${PLINK2} \
    --bfile "${TEST}" \
    --threads ${THREADS} \
    --extract "${OUT_DIR}/factor_snps_${P}.txt" \
    --score "${T2D_GWAS}" 3 6 9 header cols=+scoresums \
    --out "${OUT_DIR}/prs_factor_t2dweight_${P}"
    # VERIFY: column numbers for ID, A1 and BETA in result_t2d6.txt

done

# PLACO-based score, T2D weights
${PLINK2} \
  --bfile "${TEST}" \
  --threads ${THREADS} \
  --extract "${OUT_DIR}/placo_snps.txt" \
  --score "${T2D_GWAS}" 3 6 9 header cols=+scoresums \
  --out "${OUT_DIR}/prs_placo"

echo "PRS construction complete. Proceed to 08_prs/prs_evaluation.R"
