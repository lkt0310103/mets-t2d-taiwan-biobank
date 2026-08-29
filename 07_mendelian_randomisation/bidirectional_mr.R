# =============================================================================
# Bidirectional two-sample Mendelian randomisation
#
# Causal relationships are estimated for all pairwise combinations of the 12
# obesity- and metabolic syndrome-related indices (11 continuous traits and
# type 2 diabetes), in both directions.
#
#   Exposure instruments : training-set GWAS
#   Outcome effects      : testing-set GWAS
#
# The two sets are non-overlapping by design, so the analysis is a genuine
# two-sample design.
#
#   Instruments : p < 5e-8, clumped at r2 < 0.001 within 10,000 kb using the
#                 1000 Genomes Phase 3 EAS reference panel
#   Primary estimator : inverse-variance weighted (specified a priori)
#   Sensitivity : Cochran's Q, MR-Egger intercept, leave-one-out
#
# =============================================================================

library(TwoSampleMR)
library(data.table)

out_dir <- "results/mr"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sumstats_dir <- "data/sumstats"
EAS_PLINK    <- "reference/g1000_eas"    # local LD reference for clumping
PLINK_BIN    <- "plink"

traits <- c("bmi", "bf", "wc", "whr", "sbp", "dbp",
            "hdl", "ldl", "tg", "fg", "hba1c", "t2d")

binary_traits <- c("t2d")

# ---------------------------------------------------------------------------
# Helper: read and format summary statistics
# ---------------------------------------------------------------------------
load_sumstats <- function(trait, set, type = c("exposure", "outcome")) {

  type <- match.arg(type)
  dt <- fread(file.path(sumstats_dir, sprintf("result_%s%s.txt", trait, set)))

  dt <- dt[, .(SNP        = ID,
               beta       = BETA,
               se         = SE,
               effect_allele = A1,
               other_allele  = OMITTED,
               eaf        = A1_FREQ,
               pval       = P,
               samplesize = OBS_CT)]

  fmt <- format_data(as.data.frame(dt), type = type)
  fmt[[paste0(type, "")]] <- trait
  fmt
}

# ---------------------------------------------------------------------------
# Main loop: every ordered pair (exposure -> outcome)
# ---------------------------------------------------------------------------
results <- list()

for (exp_trait in traits) {

  # -------------------------------------------------------------------------
  # 1. Instrument selection from the training-set GWAS
  # -------------------------------------------------------------------------
  exp_dat <- load_sumstats(exp_trait, 6, "exposure")
  exp_dat <- subset(exp_dat, pval.exposure < 5e-8)

  if (nrow(exp_dat) == 0) next

  exp_dat <- clump_data(exp_dat,
                        clump_r2   = 0.001,
                        clump_kb   = 10000,
                        pop        = "EAS")
  # Local alternative, avoiding the remote API:
  # exp_dat <- ld_clump(dplyr::tibble(rsid = exp_dat$SNP,
  #                                   pval = exp_dat$pval.exposure),
  #                     clump_r2 = 0.001, clump_kb = 10000,
  #                     plink_bin = PLINK_BIN, bfile = EAS_PLINK)

  # Instrument strength: F = (beta / se)^2
  # Reported for transparency; no instrument was excluded on this basis.
  exp_dat$F_statistic <- (exp_dat$beta.exposure / exp_dat$se.exposure)^2
  message(sprintf("%s: %d instruments, F range %.1f-%.1f",
                  exp_trait, nrow(exp_dat),
                  min(exp_dat$F_statistic), max(exp_dat$F_statistic)))

  for (out_trait in traits) {

    if (out_trait == exp_trait) next

    # -----------------------------------------------------------------------
    # 2. Outcome effects from the testing-set GWAS
    # -----------------------------------------------------------------------
    out_all <- load_sumstats(out_trait, 4, "outcome")
    out_dat <- subset(out_all, SNP %in% exp_dat$SNP)

    if (nrow(out_dat) < 3) next

    # -----------------------------------------------------------------------
    # 3. Harmonise; palindromic SNPs with intermediate allele frequencies
    #    are excluded (action = 2 is the TwoSampleMR default)
    # -----------------------------------------------------------------------
    dat <- harmonise_data(exp_dat, out_dat, action = 2)
    dat <- subset(dat, mr_keep)

    if (nrow(dat) < 3) next

    # -----------------------------------------------------------------------
    # 4. Causal estimation and sensitivity analyses
    # -----------------------------------------------------------------------
    mr_res  <- mr(dat, method_list = c("mr_ivw",
                                       "mr_egger_regression",
                                       "mr_weighted_median",
                                       "mr_weighted_mode"))
    het     <- mr_heterogeneity(dat, method_list = "mr_ivw")
    pleio   <- mr_pleiotropy_test(dat)
    loo     <- mr_leaveoneout(dat)

    key <- paste(exp_trait, out_trait, sep = "_")
    results[[key]] <- list(mr = mr_res, het = het, pleio = pleio, loo = loo)

    fwrite(loo, file.path(out_dir, sprintf("loo_%s.csv", key)))
  }
}

# ---------------------------------------------------------------------------
# 5. Assemble the results table
#
# Binary outcome (type 2 diabetes): estimates are exponentiated to odds ratios
# per 1 SD increase in the genetically predicted exposure.
# ---------------------------------------------------------------------------
tab <- rbindlist(lapply(names(results), function(k) {
  r     <- results[[k]]
  ivw   <- subset(r$mr, method == "Inverse variance weighted")
  parts <- strsplit(k, "_", fixed = TRUE)[[1]]

  data.table(exposure = parts[1],
             outcome  = parts[2],
             nsnp     = ivw$nsnp,
             b        = ivw$b,
             se       = ivw$se,
             lci      = ivw$b - 1.96 * ivw$se,
             uci      = ivw$b + 1.96 * ivw$se,
             pval     = ivw$pval,
             Q_pval        = r$het$Q_pval,
             egger_int_p   = r$pleio$pval)
}), fill = TRUE)

tab[outcome %in% binary_traits,
    `:=`(OR = exp(b), OR_lci = exp(lci), OR_uci = exp(uci))]

fwrite(tab, file.path(out_dir, "mr_all_pairs.csv"))

message("Analyses completed: ", nrow(tab))
message("Nominally significant (p<0.05): ", sum(tab$pval < 0.05))

# Bonferroni threshold for the 22 type 2 diabetes-related analyses
message("Bonferroni threshold for T2D analyses: ", 0.05 / 22)

sessionInfo()
