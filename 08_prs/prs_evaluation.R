# =============================================================================
# Polygenic risk score evaluation
#
# Each standardised score is added to a base model and evaluated against the
# composite outcome of coexisting clinical metabolic syndrome and type 2
# diabetes (cases) versus neither condition (control participants).
#
#   Base model : sex, age, PC1-PC4, smoking status, alcohol use, betel nut
#                chewing, educational attainment and BMI
#   Metrics    : odds ratio per 1 SD of the PRS
#                change in discrimination (delta AUC), tested with DeLong's test
#                incremental variance explained (delta R2, Nagelkerke)
#                continuous net reclassification improvement (50 bootstraps)
#
# =============================================================================

library(data.table)
library(pROC)
library(fmsb)      # NagelkerkeR2

set.seed(123)   

pheno_file <- "data/phenotypes_testing.csv"
prs_dir    <- "results/prs"
out_dir    <- "results/prs"

N_BOOT <- 50

# ---------------------------------------------------------------------------
# 1. Outcome and covariates
# ---------------------------------------------------------------------------
dat <- fread(pheno_file)

# Composite outcome: 1 = clinical MetS AND type 2 diabetes
#                    0 = neither condition
#                    individuals meeting only one condition are excluded
dat <- dat[!is.na(clinical_mets) & !is.na(t2d)]
dat[, outcome := fifelse(clinical_mets == 1 & t2d == 1, 1L,
                  fifelse(clinical_mets == 0 & t2d == 0, 0L, NA_integer_))]
dat <- dat[!is.na(outcome)]

message("Cases: ", sum(dat$outcome == 1),
        " | Control participants: ", sum(dat$outcome == 0))

covars <- c("sex", "age", "PC1", "PC2", "PC3", "PC4",
            "smoking", "alcohol", "betelnut", "education", "BMI")

base_formula <- as.formula(
  paste("outcome ~", paste(covars, collapse = " + ")))

base_model <- glm(base_formula, data = dat, family = binomial())
base_pred  <- predict(base_model, type = "response")
base_roc   <- roc(dat$outcome, base_pred, quiet = TRUE)
base_r2    <- NagelkerkeR2(base_model)$R2

message(sprintf("Base model AUC: %.4f (%.4f, %.4f)",
                as.numeric(auc(base_roc)),
                as.numeric(ci.auc(base_roc))[1],
                as.numeric(ci.auc(base_roc))[3]))

# ---------------------------------------------------------------------------
# 2. Continuous net reclassification improvement
#
# cNRI = P(up | case) - P(down | case) + P(down | control) - P(up | control)
# ---------------------------------------------------------------------------
cnri <- function(y, p_old, p_new) {
  d <- p_new - p_old
  case <- y == 1
  ctrl <- y == 0
  (mean(d[case] > 0) - mean(d[case] < 0)) +
  (mean(d[ctrl] < 0) - mean(d[ctrl] > 0))
}

cnri_boot <- function(y, p_old, p_new, n_boot = N_BOOT) {
  n <- length(y)
  vals <- replicate(n_boot, {
    idx <- sample.int(n, n, replace = TRUE)
    cnri(y[idx], p_old[idx], p_new[idx])
  })
  c(est = cnri(y, p_old, p_new),
    lci = unname(quantile(vals, 0.025)),
    uci = unname(quantile(vals, 0.975)))
}

# ---------------------------------------------------------------------------
# 3. Evaluate each score
# ---------------------------------------------------------------------------
score_files <- list(
  list(label = "Factor GWAS (MetS weights)",
       files = file.path(prs_dir, paste0("prs_factor_metsweight_",
                          c("5e-8","1e-5","1e-3","0.01","0.05","0.1"), ".sscore")),
       thresholds = c("5e-8","1e-5","1e-3","0.01","0.05","0.1")),
  list(label = "Factor GWAS (T2D weights)",
       files = file.path(prs_dir, paste0("prs_factor_t2dweight_",
                          c("5e-8","1e-5","1e-3","0.01","0.05","0.1"), ".sscore")),
       thresholds = c("5e-8","1e-5","1e-3","0.01","0.05","0.1")),
  list(label = "PLACO (FDR < 0.05)",
       files = file.path(prs_dir, "prs_placo.sscore"),
       thresholds = "FDR<0.05")
)

results <- data.table()

for (grp in score_files) {
  for (k in seq_along(grp$files)) {

    sc <- fread(grp$files[k])
    setnames(sc, grep("SCORE1_AVG|SCORE1_SUM", names(sc), value = TRUE)[1], "PRS")

    d <- merge(dat, sc[, .(IID, PRS)], by = "IID")
    d[, PRS := as.numeric(scale(PRS))]        # zero mean, unit variance

    full_model <- glm(update(base_formula, . ~ . + PRS),
                      data = d, family = binomial())

    full_pred <- predict(full_model, type = "response")
    full_roc  <- roc(d$outcome, full_pred, quiet = TRUE)

    # Refit the base model on the same rows, so the two ROC curves are paired
    base_k     <- glm(base_formula, data = d, family = binomial())
    base_pred_k <- predict(base_k, type = "response")
    base_roc_k  <- roc(d$outcome, base_pred_k, quiet = TRUE)

    delong <- roc.test(base_roc_k, full_roc, method = "delong")

    or   <- exp(coef(summary(full_model))["PRS", "Estimate"])
    ci   <- exp(confint.default(full_model)["PRS", ])
    r2d  <- NagelkerkeR2(full_model)$R2 - NagelkerkeR2(base_k)$R2
    nri  <- cnri_boot(d$outcome, base_pred_k, full_pred)

    results <- rbind(results, data.table(
      model     = grp$label,
      threshold = grp$thresholds[k],
      nSNP      = NA_integer_,   # filled from the clumping logs
      OR        = or,
      OR_lci    = ci[1],
      OR_uci    = ci[2],
      AUC       = as.numeric(auc(full_roc)),
      AUC_lci   = as.numeric(ci.auc(full_roc))[1],
      AUC_uci   = as.numeric(ci.auc(full_roc))[3],
      dAUC      = as.numeric(auc(full_roc)) - as.numeric(auc(base_roc_k)),
      dAUC_p    = delong$p.value,
      dR2_pct   = 100 * r2d,
      cNRI      = nri["est"],
      cNRI_lci  = nri["lci"],
      cNRI_uci  = nri["uci"]))
  }
}

fwrite(results, file.path(out_dir, "prs_performance.csv"))
print(results)

# ---------------------------------------------------------------------------
# 4. Risk stratification by PRS decile
# ---------------------------------------------------------------------------
decile_or <- function(score_file) {
  sc <- fread(score_file)
  setnames(sc, grep("SCORE1_AVG|SCORE1_SUM", names(sc), value = TRUE)[1], "PRS")
  d <- merge(dat, sc[, .(IID, PRS)], by = "IID")
  d[, decile := cut(PRS, breaks = quantile(PRS, probs = seq(0, 1, 0.1)),
                    include.lowest = TRUE, labels = 1:10)]
  d[, decile := relevel(factor(decile), ref = "1")]

  m <- glm(update(base_formula, . ~ . + decile), data = d, family = binomial())
  cf <- coef(summary(m))
  rows <- grep("^decile", rownames(cf))
  data.table(decile = sub("decile", "", rownames(cf)[rows]),
             OR     = exp(cf[rows, "Estimate"]),
             lci    = exp(cf[rows, "Estimate"] - 1.96 * cf[rows, "Std. Error"]),
             uci    = exp(cf[rows, "Estimate"] + 1.96 * cf[rows, "Std. Error"]))
}

fwrite(decile_or(file.path(prs_dir, "prs_factor_metsweight_0.1.sscore")),
       file.path(out_dir, "decile_factor_mets.csv"))
fwrite(decile_or(file.path(prs_dir, "prs_placo.sscore")),
       file.path(out_dir, "decile_placo.csv"))

sessionInfo()
