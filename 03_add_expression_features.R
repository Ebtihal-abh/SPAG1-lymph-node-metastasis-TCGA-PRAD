# ============================================
# 03_add_expression_features.R
# Purpose:
# - Load grouped clinical dataset
# - Load Xena expression matrix (.gz), convert to numeric matrix
# - Save expression as RDS
# - Add SPAG1 expression
# - Save final feature dataset as RDS
# ============================================

library(dplyr)

# Expression data: TCGA-PRAD HiSeq V2 RNA-seq downloaded from UCSC Xena
# (https://xenabrowser.net/datapages/). Place file in data_raw/ before running.

# ---- Paths 
path_dat_grouped <- "data_processed/prad_master_grouped.rds"
path_expr_raw    <- "data_raw/TCGA.PRAD.sampleMap_HiSeqV2.gz"
path_expr_rds    <- "data_processed/prad_expression_matrix.rds"
path_out_rds     <- "data_processed/prad_master_features.rds"



dir.create("data_processed", showWarnings = FALSE)

# ---- 1) Load grouped clinical dataset ----
dat <- readRDS(path_dat_grouped)

stopifnot("patient_id" %in% names(dat))
stopifnot(nrow(dat) > 0)

# ---- 2) Load expression raw (gz) and convert to matrix ----
# Xena sampleMap format typically:
expr_raw <- read.delim(path_expr_raw, check.names = FALSE)

# Basic sanity
stopifnot(ncol(expr_raw) > 10)
stopifnot(nrow(expr_raw) > 1000)

gene_col <- expr_raw[[1]]
mat_use  <- as.matrix(expr_raw[, -1])

rownames(mat_use) <- gene_col
rm(expr_raw)

# Ensure numeric matrix
mode(mat_use) <- "numeric"

# ---- 3) Confirm IDs match clinical patient_id ----

col_ids <- colnames(mat_use)

# Create patient-level IDs from expression columns
col_patient <- substr(col_ids, 1, 12)

# Map from patient_id -> first matching column in expression
idx <- match(dat$patient_id, col_patient)

# Some may be NA if expression has different naming or missing samples
cat("Patients in dat:", nrow(dat), "\n")
cat("Matched expression columns:", sum(!is.na(idx)), "\n")
cat("Unmatched patients:", sum(is.na(idx)), "\n\n")

dropped <- dat$patient_id[is.na(idx)]

writeLines(dropped, "data_processed/unmatched_patients.txt")

# Keep only matched rows 
dat_use <- dat[!is.na(idx), , drop = FALSE]
mat_use2 <- mat_use[, idx[!is.na(idx)], drop = FALSE]

# Now columns align 1:1 with dat_use rows
stopifnot(nrow(dat_use) == ncol(mat_use2))
stopifnot(all(substr(colnames(mat_use2), 1, 12) == dat_use$patient_id))

# Save expression matrix as RDS for efficient loading in downstream scripts.
saveRDS(mat_use2, path_expr_rds)


# ---- 4) Add SPAG1 ----
dat_use$SPAG1 <- as.numeric(mat_use2["SPAG1", ])


# Quick sanity
cat("SPAG1 summary:\n")
print(summary(dat_use$SPAG1))
cat("\n")

# ---- 6) Save final dataset with features ----
saveRDS(dat_use, path_out_rds)

cat("Saved:\n")
cat("-", path_expr_rds, "\n")
cat("-", path_out_rds, "\n")
