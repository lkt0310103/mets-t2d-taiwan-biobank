# =============================================================================
# Exploratory factor analysis of the genetic correlation matrix
#
# EFA is applied to the LDSC-derived genetic correlation matrix using minimum
# residual (minres) extraction and oblique (Oblimin) rotation. Three traits
# (body fat percentage, waist-to-hip ratio and LDL-cholesterol) are dropped to
# reduce redundancy, leaving eight representative traits. The number of factors
# is decided from a scree plot with an eigenvalue >= 1 criterion.
#
# Input : results/ldsc/genetic_correlation_matrix.csv
#         (11 x 11 matrix of rg values assembled from the ldsc.py --rg logs)
# Output: results/factor/efa_loadings.csv, results/factor/scree_plot.pdf
#
# =============================================================================

library(psych)

set.seed(123)   

in_file   <- "results/ldsc/genetic_correlation_matrix.csv"
out_dir   <- "results/factor"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 1. Read the genetic correlation matrix
# ---------------------------------------------------------------------------
rg_all <- as.matrix(read.csv(in_file, row.names = 1, check.names = FALSE))

# Enforce symmetry and a unit diagonal; LDSC point estimates can fall slightly
# outside [-1, 1] for near-perfectly correlated traits.
rg_all <- (rg_all + t(rg_all)) / 2
diag(rg_all) <- 1
rg_all[rg_all >  1] <-  1
rg_all[rg_all < -1] <- -1

# ---------------------------------------------------------------------------
# 2. Retain the eight representative traits
#
# BF and WHR are dropped because they are near-perfectly or strongly
# correlated with BMI and WC; LDL-C is dropped because it showed no
# significant genetic correlation with the other traits.
# ---------------------------------------------------------------------------
keep <- c("BMI", "WC", "SBP", "DBP", "HDL-C", "TG", "FG", "HbA1c")
rg   <- rg_all[keep, keep]

# ---------------------------------------------------------------------------
# 3. Scree plot and eigenvalues
# ---------------------------------------------------------------------------
eigenvalues <- eigen(rg)$values
print(data.frame(component = seq_along(eigenvalues), eigenvalue = eigenvalues))

pdf(file.path(out_dir, "scree_plot.pdf"), width = 6, height = 5)
scree(rg, factors = TRUE, pc = TRUE)
dev.off()

n_factors <- sum(eigenvalues >= 1)
message("Factors retained (eigenvalue >= 1): ", n_factors)

# ---------------------------------------------------------------------------
# 4. Exploratory factor analysis
#
# n.obs is required so that psych can compute standard errors from a
# correlation matrix. VERIFY: the value supplied in the original analysis
# (e.g. the mean effective sample size across the eight training-set GWAS).
# ---------------------------------------------------------------------------
N_OBS <- 63000   # VERIFY

efa_fit <- fa(r        = rg,
              nfactors = n_factors,
              rotate   = "oblimin",
              fm       = "minres",
              n.obs    = N_OBS)

print(efa_fit, cut = 0.40, sort = FALSE, digits = 2)

loadings <- as.data.frame(unclass(efa_fit$loadings))
loadings$communality <- efa_fit$communality
loadings$uniqueness  <- efa_fit$uniquenesses

write.csv(loadings, file.path(out_dir, "efa_loadings.csv"))

# Proportion and cumulative variance explained
print(efa_fit$Vaccounted)

# Inter-factor correlations (oblique rotation)
print(efa_fit$Phi)

sessionInfo()
