# =============================================================================
# Gene-level pleiotropy analysis (PLACO+)
#
# MAGMA gene-level p values for each trait are converted to Z scores, and
# PLACO+ is applied to the same 55 trait pairs as in the SNP-level analysis.
#
#   Candidate pleiotropic genes        : FDR < 0.05
#   High-confidence pleiotropic genes  : p_placo < 0.05/18,469
#
# =============================================================================

library(data.table)
require(devtools)
source_url("https://github.com/RayDebashree/PLACO/blob/master/PLACO_v0.2.0.R?raw=TRUE")
source("PLACO_v0.2.0.R") 

magma_dir <- "results/magma"
out_dir   <- "results/pleiotropy/gene_level"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

traits <- c("bmi", "bf", "wc", "whr", "sbp", "dbp",
            "hdl", "ldl", "tg", "fg", "hba1c")

N_GENES <- 18469
BONF    <- 0.05 / N_GENES

# ---------------------------------------------------------------------------
# 1. Read MAGMA output and convert gene p values to Z scores
# ---------------------------------------------------------------------------
read_genes <- function(trait) {
  dt <- fread(file.path(magma_dir, paste0(trait, ".genes.out")))
  dt <- dt[, .(GENE, P)]
  # MAGMA reports a one-sided p value; ZSTAT is also available in the output
  # and is used directly where present.
  dt[, Z := qnorm(P / 2, lower.tail = FALSE)]
  setnames(dt, c("Z", "P"), paste0(c("Z_", "P_"), trait))
  dt
}

dat <- Reduce(function(x, y) merge(x, y, by = "GENE"),
              lapply(traits, read_genes))

Z <- as.matrix(dat[, paste0("Z_", traits), with = FALSE])
P <- as.matrix(dat[, paste0("P_", traits), with = FALSE])
rownames(Z) <- rownames(P) <- dat$GENE
colnames(Z) <- colnames(P) <- traits

message("Genes analysed: ", nrow(Z))

# ---------------------------------------------------------------------------
# 2. Pairwise PLACO+
# ---------------------------------------------------------------------------
all_pairs <- combn(traits, 2, simplify = FALSE)

summary_counts <- data.frame()

for (pr in all_pairs) {

  t1 <- pr[1]; t2 <- pr[2]
  message("Gene-level PLACO+: ", t1, " - ", t2)

  Zsub <- Z[, c(t1, t2)]
  Psub <- P[, c(t1, t2)]

  VarZ <- var.placo(Zsub, Psub, p.threshold = 1e-4)
  out  <- placo.plus(Z = Zsub, VarZ = VarZ)

  res <- data.table(GENE    = rownames(Zsub),
                    T.placo = out$T.placo,
                    p.placo = out$p.placo)
  res[, fdr := p.adjust(p.placo, method = "BH")]

  fwrite(res, file.path(out_dir, sprintf("placo_gene_%s_%s.txt.gz", t1, t2)),
         sep = "\t", compress = "gzip")

  summary_counts <- rbind(summary_counts,
    data.frame(trait1 = t1, trait2 = t2,
               n_fdr  = sum(res$fdr < 0.05),
               n_bonf = sum(res$p.placo < BONF)))
}

fwrite(summary_counts, file.path(out_dir, "pairwise_gene_counts.csv"))

# ---------------------------------------------------------------------------
# 3. Unique pleiotropic genes
# ---------------------------------------------------------------------------
genes_fdr  <- character(0)
genes_bonf <- character(0)

for (pr in all_pairs) {
  r <- fread(file.path(out_dir, sprintf("placo_gene_%s_%s.txt.gz", pr[1], pr[2])))
  genes_fdr  <- union(genes_fdr,  r$GENE[r$fdr < 0.05])
  genes_bonf <- union(genes_bonf, r$GENE[r$p.placo < BONF])
}

message("Unique pleiotropic genes, FDR < 0.05     : ", length(genes_fdr))
message("Unique pleiotropic genes, Bonferroni     : ", length(genes_bonf))

writeLines(genes_fdr,  file.path(out_dir, "pleiotropic_genes_fdr05.txt"))
writeLines(genes_bonf, file.path(out_dir, "pleiotropic_genes_bonferroni.txt"))

# The FDR-based gene list was submitted to the GENE2FUNC module of FUMA
# (https://fuma.ctglab.nl/) for functional enrichment; see README for the
# settings used, as FUMA is a web application and has no accompanying script.

sessionInfo()
