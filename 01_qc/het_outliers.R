# =============================================================================
# Identify heterozygosity outliers
#
# Flags individuals whose observed heterozygosity rate deviates by more than
# 3 SD from the sample mean, and writes a PLINK --remove exclusion list.
#
# Usage: Rscript het_outliers.R <input.het> <output.txt>
#   <input.het>  : output of `plink --het`
#   <output.txt> : two-column (FID IID) exclusion list
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript het_outliers.R <input.het> <output.txt>")
}

het_file <- args[1]
out_file <- args[2]

het <- read.table(het_file, header = TRUE, stringsAsFactors = FALSE)

# Observed heterozygosity rate = (N(NM) - O(HOM)) / N(NM)
het$het_rate <- (het$N.NM. - het$O.HOM.) / het$N.NM.

mu    <- mean(het$het_rate, na.rm = TRUE)
sigma <- sd(het$het_rate, na.rm = TRUE)

fail <- het[abs(het$het_rate - mu) > 3 * sigma, c("FID", "IID")]

write.table(fail, out_file,
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

message(sprintf("Mean heterozygosity rate: %.5f (SD %.5f)", mu, sigma))
message(sprintf("Individuals excluded (>3 SD): %d", nrow(fail)))
