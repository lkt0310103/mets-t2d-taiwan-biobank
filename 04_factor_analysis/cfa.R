# =============================================================================
# Second-order confirmatory factor analysis
#
# The four-domain structure identified by EFA is tested as a second-order
# model, in which a single latent metabolic syndrome (MetS) factor accounts
# for the covariance among four first-order factors: obesity, blood pressure,
# glycaemic and lipid.
#
# The model is fitted exclusively in the testing set (40% of the cohort) so
# that validation is independent of the discovery sample. Observed variables
# are standardised (Z score) before fitting, and parameters are estimated by
# maximum likelihood with robust standard errors (MLR).
#
# Input : data/phenotypes_testing.csv  (one row per participant; individual-level
#         Taiwan Biobank data, not redistributable)
# Output: results/factor/cfa_parameters.csv, results/factor/cfa_fit.csv
#
# =============================================================================

library(lavaan)

in_file <- "data/phenotypes_testing.csv"
out_dir <- "results/factor"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Prepare the observed variables
# ---------------------------------------------------------------------------
dat <- read.csv(in_file, stringsAsFactors = FALSE)

vars <- c("BMI", "WC", "SBP", "DBP", "FG", "HbA1c", "TG", "HDL_C")

# Standardise each observed variable
dat[vars] <- scale(dat[vars])

# ---------------------------------------------------------------------------
# 2. Specify the second-order model
# ---------------------------------------------------------------------------
model <- '
  # First-order factors
  Obesity  =~ BMI + WC
  BP       =~ SBP + DBP
  Glycemia =~ FG  + HbA1c
  Lipids   =~ TG  + HDL_C

  # Second-order factor
  MetS =~ Obesity + BP + Glycemia + Lipids
'

# ---------------------------------------------------------------------------
# 3. Fit
#
# missing = "listwise"  -> complete cases only
# missing = "ml"        -> full information maximum likelihood
# VERIFY which was used.
# ---------------------------------------------------------------------------
fit <- cfa(model,
           data      = dat,
           estimator = "MLR",
           std.lv    = FALSE,
           missing   = "listwise")

summary(fit, fit.measures = TRUE, standardized = TRUE)

# ---------------------------------------------------------------------------
# 4. Export parameter estimates and fit indices
# ---------------------------------------------------------------------------
params <- parameterEstimates(fit, standardized = TRUE)
write.csv(params, file.path(out_dir, "cfa_parameters.csv"), row.names = FALSE)

fit_indices <- fitMeasures(fit, c("cfi.robust", "tli.robust",
                                  "rmsea.robust",
                                  "rmsea.ci.lower.robust",
                                  "rmsea.ci.upper.robust",
                                  "srmr"))
print(fit_indices)
write.csv(as.data.frame(t(fit_indices)),
          file.path(out_dir, "cfa_fit.csv"), row.names = FALSE)

# Acceptance criteria used in the manuscript:
#   robust CFI and TLI >= 0.95, robust RMSEA < 0.08, SRMR < 0.08

# Standardised second-order loadings, carried forward as fixed values in the
# factor-based GWAS (see 05_factor_gwas/)
second_order <- subset(params, lhs == "MetS" & op == "=~",
                       select = c("rhs", "est", "se", "z", "pvalue", "std.all"))
print(second_order)

sessionInfo()
