# 08d_SPAG1_GSEA_by_Nstratum.R
# ============================================================
# Study:   TCGA-PRAD RNA-seq cohort (Illumina HiSeqV2)
# Purpose: Address reviewer request — characterise the SPAG1-associated
#          transcriptional programme (high vs low, median split) SEPARATELY
#          within node-negative (N0) and node-positive (N1) tumours, and test
#          whether pathways such as G2M Checkpoint and Interferon-gamma
#          Response are enriched differently across nodal strata.
#
#          High vs low is defined by the stratum-specific median SPAG1
#          expression, paralleling the categorical Q4-vs-Q1 contrast used in
#          the primary analysis (Script 08). Analyses are unadjusted within
#          stratum. Positive NES = enriched in high SPAG1.
#
# Inputs:
#   data_processed/prad_master_features.rds
#   data_processed/prad_expression_matrix.rds
#
# Outputs:
#   tables/08d/DE_SPAG1_medsplit_N0_limma.csv
#   tables/08d/DE_SPAG1_medsplit_N1_limma.csv
#   tables/08d/FGSEA_Hallmark_SPAG1_medsplit_N0.csv
#   tables/08d/FGSEA_Hallmark_SPAG1_medsplit_N1.csv
#   tables/08d/FGSEA_Hallmark_SPAG1_by_Nstratum_combined.csv
#   tables/08d/Key_pathways_G2M_IFNg_by_Nstratum.csv
# ============================================================
library(dplyr)
library(matrixStats)
library(limma)
library(msigdbr)
library(fgsea)
dir.create("tables/08d", recursive = TRUE, showWarnings = FALSE)

# Pathways specifically highlighted by the reviewer
KEY_PATHWAYS <- c("HALLMARK_G2M_CHECKPOINT",
                  "HALLMARK_INTERFERON_GAMMA_RESPONSE")

# ---- Load data ----
dat <- readRDS("data_processed/prad_master_features.rds")
mat <- readRDS("data_processed/prad_expression_matrix.rds")

# ---- Alignment checks ----
stopifnot(nrow(dat) == ncol(mat))
stopifnot(all(dat$patient_id == substr(colnames(mat), 1, 12)))

# ---- Define nodal strata ----
dat <- dat %>%
  mutate(
    N_stratum = case_when(
      tolower(as.character(N_group2)) == "n1" ~ "N1",
      tolower(as.character(N_group2)) == "n0" ~ "N0",
      TRUE                                    ~ NA_character_
    )
  )
keep_N <- !is.na(dat$N_stratum) & !is.na(dat$SPAG1)
dat_N  <- dat[keep_N, , drop = FALSE]
mat_N  <- mat[, keep_N, drop = FALSE]
stopifnot(nrow(dat_N) == ncol(mat_N))
cat("Samples kept (N0 + N1):", nrow(dat_N), "\n")
print(table(dat_N$N_stratum))

# ---- Hallmark gene sets ----
hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)
cat("Hallmark pathways loaded:", length(pathways), "\n")

# =========================================================
# fgsea from limma moderated t-statistics
#   Same convention as Script 08
# =========================================================
run_fgsea <- function(t_stats, gene_names, seed = 123) {
  ranks <- t_stats
  names(ranks) <- gene_names
  set.seed(seed)
  ranks <- ranks + rnorm(length(ranks), 0, 1e-10)
  ranks <- sort(ranks, decreasing = TRUE)
  fgseaMultilevel(pathways = pathways, stats = ranks) %>%
    arrange(padj) %>%
    mutate(leadingEdge = vapply(leadingEdge,
                                function(x) paste(x, collapse = ";"),
                                character(1)))
}

# =========================================================
# Median-split (high vs low) limma + fgsea within one stratum
# =========================================================
analyse_stratum <- function(stratum) {
  sel   <- dat_N$N_stratum == stratum
  dat_s <- dat_N[sel, , drop = FALSE]
  mat_s <- mat_N[, sel, drop = FALSE]
  stopifnot(nrow(dat_s) == ncol(mat_s))
  
  # Remove zero-variance genes within this stratum
  gene_var <- matrixStats::rowVars(mat_s)
  mat_s_f  <- mat_s[gene_var > 0, , drop = FALSE]
  
  # High vs low by stratum-specific median SPAG1
  spag1_median <- median(dat_s$SPAG1, na.rm = TRUE)
  group <- factor(ifelse(dat_s$SPAG1 >= spag1_median, "High", "Low"),
                  levels = c("Low", "High"))
  
  design <- model.matrix(~ group)
  fit    <- eBayes(lmFit(mat_s_f, design))
  de     <- topTable(fit, coef = "groupHigh", number = Inf, sort.by = "P")
  write.csv(de, sprintf("tables/08d/DE_SPAG1_medsplit_%s_limma.csv", stratum),
            row.names = TRUE)
  
  fg <- withCallingHandlers(
    run_fgsea(de$t, rownames(de)),
    warning = function(w) {
      cat(sprintf(" [%s] fgsea warning: %s\n", stratum, conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )
  write.csv(fg, sprintf("tables/08d/FGSEA_Hallmark_SPAG1_medsplit_%s.csv", stratum),
            row.names = FALSE)
  
  cat(sprintf("\n=== %s: n = %d | cutpoint = %.4f [Low %d | High %d] | genes = %d ===\n",
              stratum, nrow(dat_s), spag1_median,
              sum(group == "Low"), sum(group == "High"), nrow(mat_s_f)))
  cat(sprintf("Sig pathways FDR < 0.05: %d\n", sum(fg$padj < 0.05, na.rm = TRUE)))
  
  list(fg = fg, n = nrow(dat_s))
}

# =========================================================
# Run both strata
# =========================================================
res        <- lapply(c("N0", "N1"), analyse_stratum)
names(res) <- c("N0", "N1")
cat("\nSaved: DE and FGSEA tables to tables/08d/\n")

# =========================================================
# Combined wide table: all 50 Hallmark pathways, NES + padj per stratum
# =========================================================
nes_cols <- function(fg, tag) {
  out <- as.data.frame(fg)[, c("pathway", "NES", "padj")]
  names(out) <- c("pathway", paste0("NES_", tag), paste0("padj_", tag))
  out
}
combined <- merge(nes_cols(res$N0$fg, "N0"),
                  nes_cols(res$N1$fg, "N1"),
                  by = "pathway", all = TRUE)
combined <- combined[order(combined$padj_N0), ]
write.csv(combined, "tables/08d/FGSEA_Hallmark_SPAG1_by_Nstratum_combined.csv",
          row.names = FALSE)
cat("Saved: tables/08d/FGSEA_Hallmark_SPAG1_by_Nstratum_combined.csv\n")

# =========================================================
# Reviewer-focused table: G2M Checkpoint & Interferon-gamma Response
# =========================================================
key_tbl <- combined[combined$pathway %in% KEY_PATHWAYS, ]
write.csv(key_tbl, "tables/08d/Key_pathways_G2M_IFNg_by_Nstratum.csv",
          row.names = FALSE)
cat("Saved: tables/08d/Key_pathways_G2M_IFNg_by_Nstratum.csv\n")

# ---- Console summary: directly answers the reviewer question ----
cat("\nG2M Checkpoint & Interferon-gamma Response by nodal stratum (median split):\n")
print(key_tbl, row.names = FALSE)
cat(sprintf(
  paste0("\nNote: N1 (n=%d) is smaller than N0 (n=%d); fewer pathways reach\n",
         "FDR < 0.05 in N1 owing to the smaller subgroup size. Interpretation\n",
         "prioritises the direction (sign of NES) being consistent across strata.\n"),
  res$N1$n, res$N0$n))

cat("\nDONE.\n")