# Shared constants and helper functions used across multiple scripts

# TSS collapsing threshold — sites with fewer than this many samples
# are collapsed into "Other". Used in Scripts 08b, 08c and 13.
TSS_MIN_N <- 25

collapse_tss <- function(tss_vec, min_n = TSS_MIN_N) {
  counts <- table(tss_vec)
  keep   <- names(counts[counts >= min_n])
  droplevels(factor(ifelse(tss_vec %in% keep, as.character(tss_vec), "Other")))
}