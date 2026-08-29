# =============================================================================
# SNP-level pleiotropy analysis (PLACO+)
#
# Pleiotropy is tested for each of the 55 pairwise combinations of the 11
# continuous obesity- and metabolic syndrome-related traits, using the product
# of Z statistics under a composite null hypothesis. PLACO+ extends the
# original PLACO method through an inflated variance model that accommodates
# variants associated with neither trait, or with only one trait, under the
# null.
#
# Primary significance threshold : FDR < 0.05 (Benjamini-Hochberg)
# Comparison threshold           : p_placo < 5e-8
#
# The FDR-based SNP set is carried forward to functional annotation and to the
# PLACO-based polygenic risk score.
#
# =============================================================================

library(data.table)

require(devtools)
source_url("https://github.com/RayDebashree/PLACO/blob/master/PLACO_v0.2.0.R?raw=TRUE")
source("PLACO_v0.2.0.R") 
sumstats_dir <- "data/sumstats"
out_dir      <- "results/pleiotropy/snp_level"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

traits <- c("bmi", "bf", "wc", "whr", "sbp", "dbp",
            "hdl", "ldl", "tg", "fg", "hba1c")

# ---------------------------------------------------------------------------
# 1. Read training-set summary statistics and build Z / P matrices
# ---------------------------------------------------------------------------
read_trait <- function(trait) {
  dt <- fread(file.path(sumstats_dir, paste0("result_", trait, "6.txt")))
  # VERIFY: column names of your PLINK 2.0 --glm output
  dt <- dt[, .(SNP = ID, Z = BETA / SE, P = P)]
  setnames(dt, c("Z", "P"), paste0(c("Z_", "P_"), trait))
  dt
}

dat <- Reduce(function(x, y) merge(x, y, by = "SNP"),
              lapply(traits, read_trait))

Z <- as.matrix(dat[, paste0("Z_", traits), with = FALSE])
P <- as.matrix(dat[, paste0("P_", traits), with = FALSE])
rownames(Z) <- rownames(P) <- dat$SNP
colnames(Z) <- colnames(P) <- traits

# ---------------------------------------------------------------------------
# 2. Decorrelate Z statistics
#
# The traits are analysed in the same sample and are correlated, so the Z
# statistics are decorrelated before testing.
# ---------------------------------------------------------------------------
R <- cor.pearson(Z, P, p.threshold = 1e-4)   # correlation from null variants
Z_dec <- decorrelate.Z(Z, R)

# ---------------------------------------------------------------------------
# 3. Pairwise PLACO+
# ---------------------------------------------------------------------------
all_pairs <- combn(traits, 2, simplify = FALSE)
stopifnot(length(all_pairs) == 55)

summary_counts <- data.frame()

for (pr in all_pairs) {

  t1 <- pr[1]; t2 <- pr[2]
  message("PLACO+: ", t1, " - ", t2)

  Zsub <- Z_dec[, c(t1, t2)]
  Psub <- P[,     c(t1, t2)]

  # Variance parameters estimated under the composite null
  VarZ <- var.placo(Zsub, Psub, p.threshold = 1e-4)

  out <- placo.plus(Z = Zsub, VarZ = VarZ)

  res <- data.table(SNP    = rownames(Zsub),
                    Z1     = Zsub[, 1],
                    Z2     = Zsub[, 2],
                    T.placo = out$T.placo,
                    p.placo = out$p.placo)

  res[, fdr := p.adjust(p.placo, method = "BH")]

  fwrite(res, file.path(out_dir, sprintf("placo_%s_%s.txt.gz", t1, t2)),
         sep = "\t", compress = "gzip")

  summary_counts <- rbind(summary_counts,
    data.frame(trait1 = t1, trait2 = t2,
               n_fdr        = sum(res$fdr < 0.05),
               n_genomewide = sum(res$p.placo < 5e-8)))
}

fwrite(summary_counts, file.path(out_dir, "pairwise_counts.csv"))

# ---------------------------------------------------------------------------
# 4. Unique pleiotropic SNPs across the 55 pairs
#
# The FDR-based set is the primary result and the input to the PLACO-based PRS.
# ---------------------------------------------------------------------------
collect <- function(threshold_col, cutoff) {
  ids <- character(0)
  for (pr in all_pairs) {
    f <- file.path(out_dir, sprintf("placo_%s_%s.txt.gz", pr[1], pr[2]))
    r <- fread(f)
    ids <- union(ids, r$SNP[r[[threshold_col]] < cutoff])
  }
  ids
}

snps_fdr        <- collect("fdr", 0.05)
snps_genomewide <- collect("p.placo", 5e-8)

message("Unique pleiotropic SNPs, FDR < 0.05        : ", length(snps_fdr))
message("Unique pleiotropic SNPs, p < 5e-8          : ", length(snps_genomewide))

writeLines(snps_fdr,        file.path(out_dir, "pleiotropic_snps_fdr05.txt"))
writeLines(snps_genomewide, file.path(out_dir, "pleiotropic_snps_gw.txt"))

sessionInfo()
