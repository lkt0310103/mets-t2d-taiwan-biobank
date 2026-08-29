# =============================================================================
# Factor-based genome-wide association study of the latent MetS factor
#
# The validated second-order model is re-specified as a GWAS model by adding a
# structural path from each SNP to the second-order MetS factor. Factor
# loadings are held fixed at their previously estimated values so that only the
# SNP-MetS path is re-estimated at each variant, and a lower-bound constraint
# is imposed on the residual variance of BMI to prevent an improper (Heywood)
# solution.
#
# Software: GenomicSEM (https://github.com/GenomicSEM/GenomicSEM)
# Reference: 1000 Genomes Phase 3 EAS; HapMap3 SNP list
#
# =============================================================================

library(GenomicSEM)

out_dir <- "results/factor_gwas"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

REF_EAS   <- "reference/eas_ldscores/"
HM3       <- "reference/w_hm3.snplist"
REF_PANEL <- "reference/1000G_EAS.txt"   # reference file for sumstats()

traits <- c("BMI", "WC", "SBP", "DBP", "HDL_C", "TG", "FG", "HbA1c")
files  <- file.path("data/sumstats",
                    paste0("result_", tolower(traits), "6.txt"))


sample_sizes <- c(63081, 63287, 63242, 63354, 63306, 62967, 62601, 62569)

# ---------------------------------------------------------------------------
# 1. Munge summary statistics to the HapMap3 SNP list
# ---------------------------------------------------------------------------
munge(files       = files,
      hm3         = HM3,
      trait.names = traits,
      N           = sample_sizes,
      info.filter = 0.8,
      maf.filter  = 0.01)

# ---------------------------------------------------------------------------
# 2. Multivariable LDSC -- genetic and sampling covariance matrices
# ---------------------------------------------------------------------------
LDSCoutput <- ldsc(traits      = paste0(traits, ".sumstats.gz"),
                   sample.prev = rep(NA, length(traits)),   # continuous traits
                   population.prev = rep(NA, length(traits)),
                   ld          = REF_EAS,
                   wld         = REF_EAS,
                   trait.names = traits)

saveRDS(LDSCoutput, file.path(out_dir, "ldsc_covstruct.rds"))

# ---------------------------------------------------------------------------
# 3. Prepare SNP-level summary statistics
#
# Standardised SNP effect sizes are obtained assuming an additive linear model
# for all traits (linprob = FALSE, OLS = TRUE for continuous traits).
# ---------------------------------------------------------------------------
sumstats_out <- sumstats(files        = files,
                         ref          = REF_PANEL,
                         trait.names  = traits,
                         se.logit     = rep(FALSE, length(traits)),
                         OLS          = rep(TRUE,  length(traits)),
                         linprob      = rep(FALSE, length(traits)),
                         N            = sample_sizes,
                         info.filter  = 0.8,
                         maf.filter   = 0.01)

saveRDS(sumstats_out, file.path(out_dir, "sumstats_prepared.rds"))

# ---------------------------------------------------------------------------
# 4. Specify the GWAS model
#
# First- and second-order loadings are fixed at the values estimated in the
# testing-set CFA (04_factor_analysis/cfa.R). Only the SNP -> MetS path is
# estimated at each variant.
#
# ---------------------------------------------------------------------------
model <- '
  # First-order factors, loadings fixed at CFA estimates
  Obesity  =~ 0.883*BMI  + 0.947*WC
  BP       =~ 0.880*SBP  + 0.875*DBP
  Glycemia =~ 0.791*FG   + 0.779*HbA1c
  Lipids   =~ 0.695*TG   + -0.703*HDL_C

  # Second-order factor, loadings fixed at CFA estimates
  MetS =~ 0.778*Obesity + 0.553*BP + 0.511*Glycemia + 0.761*Lipids

  # SNP effect on the second-order factor (the parameter of interest)
  MetS ~ SNP

  # Lower-bound constraint on the residual variance of BMI, to prevent an
  # improper (negative variance) solution
  BMI ~~ a*BMI
  a > 0.001
'

# ---------------------------------------------------------------------------
# 5. Run the factor GWAS
#
# sub = c("MetS ~ SNP") restricts the output to the parameter of interest.
# Set cores to the number available; parallel = TRUE is strongly recommended.
# ---------------------------------------------------------------------------
factor_gwas <- userGWAS(covstruc = LDSCoutput,
                        SNPs     = sumstats_out,
                        model    = model,
                        sub      = c("MetS ~ SNP"),
                        estimation = "DWLS",
                        parallel = TRUE,
                        cores    = 20,
                        toler    = 1e-40)

res <- factor_gwas[[1]]
write.csv(res, file.path(out_dir, "factor_gwas_all.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 6. Genome-wide significant loci
# ---------------------------------------------------------------------------
sig <- subset(res, Pval_Estimate < 5e-8)
message("Genome-wide significant SNPs: ", nrow(sig))
write.csv(sig, file.path(out_dir, "factor_gwas_significant.csv"), row.names = FALSE)

sessionInfo()
