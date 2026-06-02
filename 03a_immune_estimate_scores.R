# ============================================
# 03a_immune_estimate.R
# Purpose:
# - Run ESTIMATE on existing expression matrix
#
# Notes:
# - This script only computes and saves ESTIMATE scores.
# - Correlation analyses
# - All 497 patients are confirmed primary tumor (Script 01 filters
#   to sample_type == "Primary Tumor" before any downstream steps)
# - TumorPurity from ESTIMATE not used: ESTIMATE purity formula is
#   calibrated for microarray data, not RNA-seq. Purity is addressed
#   separately using TIMER2
# ============================================

library(estimate)
library(dplyr)

# ── 1. Paths ──────────────────────────────────────────────────────────────────
expr_path <- "data_processed/prad_expression_matrix.rds"
out_rds   <- "data_processed/estimate_scores.rds"
out_dir   <- "tables/03a/"
master_path  <- "data_processed/prad_master_features.rds"
enriched_rds <- "data_processed/prad_master_features.rds"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── 2. Load expression matrix ─────────────────────────────────────────────────
cat("Loading expression matrix...\n")
mat <- readRDS(expr_path)
cat("Dimensions:", nrow(mat), "genes x", ncol(mat), "patients\n")
# Expected: 20530 genes x 497 patients

# ── 3. Write expression matrix to temp file for ESTIMATE ──────────────────────
tmp_expr   <- tempfile(fileext = ".txt")
tmp_genes  <- tempfile(fileext = ".txt")
tmp_scores <- paste0(out_dir, "ESTIMATE_raw_output.txt")

mat_df        <- as.data.frame(mat)
mat_df        <- cbind(GeneSymbol = rownames(mat), mat_df)

write.table(
  mat_df,
  tmp_expr,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)
cat("Expression matrix written to temp file\n")

# ── 4. Run ESTIMATE ───────────────────────────────────────────────────────────
# Step 1: filter to ESTIMATE's internal 10,412-gene list
filterCommonGenes(
  input.f  = tmp_expr,
  output.f = tmp_genes,
  id       = "GeneSymbol"
)

# Step 2: compute ImmuneScore and StromalScore
# platform = "illumina" is correct for TCGA HiSeq RNA-seq data
estimateScore(
  input.ds  = tmp_genes,
  output.ds = tmp_scores,
  platform  = "illumina"
)
cat("ESTIMATE scoring complete\n")

# ── 5. Read and clean ESTIMATE output ─────────────────────────────────────────
# ESTIMATE writes GCT format: 2 header lines to skip
# check.names = FALSE prevents R converting dashes to dots in column names
raw <- read.table(
  tmp_scores,
  header      = TRUE,
  sep         = "\t",
  skip        = 2,
  row.names   = 1,
  check.names = FALSE
)

# Remove Description column if present
raw <- raw[, colnames(raw) != "Description", drop = FALSE]
cat("ESTIMATE output:", nrow(raw), "score rows x", ncol(raw), "patients\n")
# Expected rows: StromalScore, ImmuneScore, ESTIMATEScore

# ── 6. Transpose and clean patient IDs ───────────────────────────────────────
scores_df             <- as.data.frame(t(raw))
scores_df$barcode_raw <- rownames(scores_df)

# Safety: convert any dots back to dashes
# (R sometimes converts dashes to dots in rownames even with check.names=FALSE)
scores_df$barcode_raw <- gsub("\\.", "-", scores_df$barcode_raw)

# Trim to 12-character patient ID — matches format used throughout pipeline
scores_df$patient_id <- substr(scores_df$barcode_raw, 1, 12)

# Check for duplicates after trimming
n_dupes <- sum(duplicated(scores_df$patient_id))
if (n_dupes > 0) {
  cat("Note:", n_dupes, "duplicate patient IDs after trimming — keeping first\n")
  scores_df <- scores_df[!duplicated(scores_df$patient_id), ]
}
cat("ESTIMATE scores ready for", nrow(scores_df), "patients\n")
# Expected: 497

# ── 7. Save scores ────────────────────────────────────────────────────────────
score_cols <- intersect(
  c("StromalScore", "ImmuneScore", "ESTIMATEScore"),
  colnames(scores_df)
)
cat("Score columns:", paste(score_cols, collapse = ", "), "\n")

scores_out <- scores_df[, c("patient_id", score_cols)]

# CSV for human inspection
write.csv(
  scores_out,
  paste0(out_dir, "ESTIMATE_patient_scores.csv"),
  row.names = FALSE
)

# Save RDS 
saveRDS(scores_out, out_rds)

cat("Scores saved to:\n")
cat(" ", paste0(out_dir, "ESTIMATE_patient_scores.csv"), "\n")
cat(" ", out_rds, "\n")

# ── 8. Merge into master dataset and save ───────────────────────────────────────
master <- readRDS(master_path)
master_immune <- master %>%
  left_join(scores_out[, c("patient_id", "ImmuneScore", "StromalScore")],
            by = "patient_id")
saveRDS(master_immune, enriched_rds)
cat("Enriched master dataset saved to:", enriched_rds, "\n")
cat("Columns added: ImmuneScore, StromalScore\n")

# ── 9. Final summary ──────────────────────────────────────────────────────────
cat("\n════ Script 03a Complete ════\n")
cat("ESTIMATE scores computed for", nrow(scores_out), "patients\n")
cat("Columns saved:", paste(score_cols, collapse = ", "), "\n")
cat("Enriched master RDS saved to:", enriched_rds, "\n")
cat("All downstream scripts load prad_master_features.rds which now includes ImmuneScore\n")


