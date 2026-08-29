# =============================================================================
# Build the exclusion list for cryptic relatedness
#
# For each pair of individuals with PI_HAT >= 0.1875 (third-degree or closer),
# the member with the higher genotype missingness rate is excluded. Pairs are
# processed iteratively so that an individual already marked for exclusion is
# not counted again.
#
# Usage: Rscript related_exclusion.R <related.genome> <miss.imiss> <output.txt>
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript related_exclusion.R <related.genome> <miss.imiss> <output.txt>")
}

genome_file <- args[1]
imiss_file  <- args[2]
out_file    <- args[3]

genome <- read.table(genome_file, header = TRUE, stringsAsFactors = FALSE)
imiss  <- read.table(imiss_file,  header = TRUE, stringsAsFactors = FALSE)

# Missingness rate keyed by FID_IID
miss_rate <- setNames(imiss$F_MISS, paste(imiss$FID, imiss$IID, sep = "_"))

id1 <- paste(genome$FID1, genome$IID1, sep = "_")
id2 <- paste(genome$FID2, genome$IID2, sep = "_")

# Process pairs in descending order of relatedness
ord <- order(genome$PI_HAT, decreasing = TRUE)

excluded <- character(0)

for (i in ord) {
  a <- id1[i]
  b <- id2[i]

  # Skip if either member has already been removed
  if (a %in% excluded || b %in% excluded) next

  ma <- ifelse(a %in% names(miss_rate), miss_rate[[a]], NA_real_)
  mb <- ifelse(b %in% names(miss_rate), miss_rate[[b]], NA_real_)

  # Exclude the member with higher missingness; ties broken by the first member
  drop <- if (is.na(ma) || is.na(mb)) a else if (mb > ma) b else a
  excluded <- c(excluded, drop)
}

if (length(excluded) == 0L) {
  fail <- data.frame(FID = character(0), IID = character(0))
} else {
  parts <- do.call(rbind, strsplit(excluded, "_", fixed = TRUE))
  fail  <- data.frame(FID = parts[, 1], IID = parts[, 2],
                      stringsAsFactors = FALSE)
}

write.table(fail, out_file,
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

message(sprintf("Related pairs (PI_HAT >= 0.1875): %d", nrow(genome)))
message(sprintf("Individuals excluded: %d", nrow(fail)))
