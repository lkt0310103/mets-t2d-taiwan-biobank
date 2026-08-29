#!/bin/bash
# =============================================================================
# Gene-based association analysis (MAGMA)
#
# SNP-level summary statistics for each trait are aggregated into gene-level
# statistics across 18,469 protein-coding genes (NCBI Build 38 annotation),
# with intragenic linkage disequilibrium modelled using the 1000 Genomes
# Phase 3 East Asian reference panel. Gene-level statistics are produced for:
#
#   - the 11 continuous traits, as input to the gene-level PLACO+ analysis
#   - the latent MetS factor GWAS, for the gene-based results in the main text
#
# Software: MAGMA v1.10   <-- VERIFY version
# =============================================================================

set -euo pipefail

MAGMA="./magma/magma"
ANNOT_DIR="./reference/magma"
REF_EAS="${ANNOT_DIR}/g1000_eas"              # 1000G EAS PLINK binary
GENE_LOC="${ANNOT_DIR}/NCBI38.gene.loc"       # protein-coding gene locations
SUMSTATS_DIR="./data/sumstats"
OUT_DIR="./results/magma"

mkdir -p "${OUT_DIR}"

TRAITS="bmi bf wc whr sbp dbp hdl ldl tg fg hba1c factor"

# ---------------------------------------------------------------------------
# 1. Annotate SNPs to genes
#
# VERIFY: whether a window around each gene was used (e.g. --annotate window=0)
# ---------------------------------------------------------------------------
${MAGMA} \
  --annotate \
  --snp-loc "${ANNOT_DIR}/snp_loc_hg38.txt" \
  --gene-loc "${GENE_LOC}" \
  --out "${OUT_DIR}/annotation"

# ---------------------------------------------------------------------------
# 2. Gene-based analysis for each trait
# ---------------------------------------------------------------------------
for TRAIT in ${TRAITS}; do

  # VERIFY: --pval column names and the sample size argument (N= or ncol=)
  ${MAGMA} \
    --bfile "${REF_EAS}" \
    --pval "${SUMSTATS_DIR}/result_${TRAIT}6.txt" use=ID,P ncol=OBS_CT \
    --gene-annot "${OUT_DIR}/annotation.genes.annot" \
    --out "${OUT_DIR}/${TRAIT}"

done

# ---------------------------------------------------------------------------
# 3. Significance thresholds
#
#   Gene-based analysis of the latent MetS factor : p < 0.05/18,469
#   Gene-level PLACO+                             : FDR < 0.05 (candidate)
#                                                   p < 0.05/18,469 (high confidence)
# ---------------------------------------------------------------------------
echo "Bonferroni threshold: $(python3 -c 'print(0.05/18469)')"

echo "MAGMA complete. Gene-level results in ${OUT_DIR}"
