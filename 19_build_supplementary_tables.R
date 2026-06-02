# =====================================================================
# 19_build_supplementary_tables.R
# ---------------------------------------------------------------------
# Submission-ready supplementary tables.
#
# Inputs: CSV files from tables/08/, tables/08b/, tables/11/, tables/11a/,
#         tables/12_validation/, tables/13/, tables/15/
#         Microenvironment correlations (ESTIMATE + TIMER2.0)
# Output: tables/Supplementary_Tables.xlsx
#         (one sheet per supplementary table)
# =====================================================================

library(openxlsx)
library(dplyr)

# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------
safe_read <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("File not found: %s", path), call. = FALSE)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}

# Format p-values: scientific notation if < 0.001, else 3 decimals
format_pvalue <- function(x) {
  ifelse(is.na(x), NA_character_,
         ifelse(x < 0.001,
                format(x, scientific = TRUE, digits = 2),
                format(round(x, 3), nsmall = 3)))
}

# Format a data frame for display:

format_table <- function(df, digits = 3) {
  for (col in names(df)) {
    if (is.numeric(df[[col]])) {
      # Detect p-value column by name pattern 
      is_pval <- grepl("p[._]val|^p$|^pval|^padj|adj.*p|fdr|wilcoxon|lrt.*p|_p_|_p$",
                       col, ignore.case = TRUE)
      
      if (is_pval) {
        df[[col]] <- format_pvalue(df[[col]])
      } else {
        df[[col]] <- round(df[[col]], digits)
      }
    }
  }
  df
}

`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

extract_or_row <- function(or, ci_low, ci_high, p, label, n_value = NA) {
  data.frame(
    Specification = label,
    n             = n_value,
    OR            = round(or, 3),
    CI_low        = round(ci_low, 3),
    CI_high       = round(ci_high, 3),
    p_value       = format(p, digits = 3, scientific = TRUE),
    stringsAsFactors = FALSE
  )
}
# ---------------------------------------------------------------------
# Workbook + styles
# ---------------------------------------------------------------------
wb <- createWorkbook()

hdr_style <- createStyle(textDecoration = "bold",
                         fgFill = "#D9E1F2",
                         halign = "center",
                         border = "Bottom")

title_style <- createStyle(textDecoration = "bold", fontSize = 12)

# =====================================================================
# Supp Table S1 — Model selection and parameterisation
# =====================================================================
addWorksheet(wb, "S1_Model_selection")

mc <- safe_read("tables/11/Logistic_N1_primary_model.csv")
mc <- format_table(mc, 3)

writeData(wb, "S1_Model_selection",
          "Supplementary Table S1. Logistic regression model selection and parameterisation",
          startRow = 1, startCol = 1)
addStyle(wb, "S1_Model_selection", title_style, rows = 1, cols = 1)

writeData(wb, "S1_Model_selection",
          "Panel A. Primary model comparison (n = 420)",
          startRow = 3, startCol = 1)
addStyle(wb, "S1_Model_selection", createStyle(textDecoration = "bold"),
         rows = 3, cols = 1)

writeData(wb, "S1_Model_selection", mc, startRow = 4, startCol = 1,
          headerStyle = hdr_style)

qr <- safe_read("tables/11/SPAG1_quartile_robustness.csv")
qr <- format_table(qr, 3)

panel_b_row <- 4 + nrow(mc) + 3
writeData(wb, "S1_Model_selection",
          "Panel B. SPAG1 parameterisation robustness",
          startRow = panel_b_row, startCol = 1)
addStyle(wb, "S1_Model_selection", createStyle(textDecoration = "bold"),
         rows = panel_b_row, cols = 1)

writeData(wb, "S1_Model_selection", qr, startRow = panel_b_row + 1,
          startCol = 1, headerStyle = hdr_style)

# =====================================================================
# Supp Table S2 — Sensitivity analyses
# =====================================================================
addWorksheet(wb, "S2_Sensitivity")

primary  <- safe_read("tables/11/SPAG1_OR_all_models.csv")
nx_n0    <- safe_read("tables/11a/Sensitivity_NX_as_N0_SPAG1.csv")
tss_stab <- safe_read("tables/13/SPAG1_stability_main_vs_TSS.csv")

primary_spag1 <- primary[
  primary$Model == "Adjusted (T stage + Gleason + SPAG1)", , drop = FALSE
]

nx_spag1 <- nx_n0[
  grepl("SPAG1", nx_n0$term, ignore.case = TRUE), , drop = FALSE
]

tss_row <- tss_stab[
  grepl("Sensitivity", tss_stab$Model), , drop = FALSE
]

sens_tbl <- rbind(
  extract_or_row(
    or      = primary_spag1$OR,
    ci_low  = primary_spag1$CI_low,
    ci_high = primary_spag1$CI_high,
    p       = as.numeric(primary_spag1$p_value),
    label   = "Primary (T + Gleason + SPAG1)",
    n_value = 420
  ),
  extract_or_row(
    or      = nx_spag1$estimate,
    ci_low  = nx_spag1$conf.low,
    ci_high = nx_spag1$conf.high,
    p       = nx_spag1$p.value,
    label   = "NX patients reclassified as N0",
    n_value = 493
  ),
  extract_or_row(
    or      = tss_row$OR,
    ci_low  = tss_row$CI_low,
    ci_high = tss_row$CI_high,
    p       = tss_row$p_value,
    label   = "Additional adjustment for tissue source site",
    n_value = 420
  )
)
writeData(wb, "S2_Sensitivity",
          "Supplementary Table S2. Robustness of the SPAG1-N1 association across pre-specified sensitivity analyses",
          startRow = 1, startCol = 1)
addStyle(wb, "S2_Sensitivity", title_style, rows = 1, cols = 1)
writeData(wb, "S2_Sensitivity", sens_tbl, startRow = 3, startCol = 1,
          headerStyle = hdr_style)
# =====================================================================
# Supp Table S3 — SPAG1 fgsea across specifications
# =====================================================================
addWorksheet(wb, "S3_SPAG1_fgsea")

spec1 <- safe_read("tables/08/FGSEA_Hallmark_SPAG1_Q4_vs_Q1_multilevel.csv")
spec1$Specification <- "SPAG1 Q4 vs Q1"
all_fgsea <- spec1
all_fgsea <- all_fgsea[, c("Specification",
                           setdiff(names(all_fgsea), "Specification"))]
all_fgsea <- format_table(all_fgsea, 4)
writeData(wb, "S3_SPAG1_fgsea",
          "Supplementary Table S3. Hallmark fgsea results — SPAG1 Q4 vs Q1",
          startRow = 1, startCol = 1)
addStyle(wb, "S3_SPAG1_fgsea", title_style, rows = 1, cols = 1)
writeData(wb, "S3_SPAG1_fgsea", all_fgsea, startRow = 3, startCol = 1,
          headerStyle = hdr_style)
# =====================================================================
# Supp Table S4 — N1 transcriptional signature and concordance
# =====================================================================
addWorksheet(wb, "S4a_N1_DE_genes")
addWorksheet(wb, "S4b_N1_fgsea")
addWorksheet(wb, "S4c_Concordance")

de_n1 <- safe_read("tables/08b/DE_N1_vs_N0_adjTGleason_TSS_limma.csv")
de_n1_sig <- de_n1[de_n1$adj.P.Val < 0.05, ]
de_n1_sig <- format_table(de_n1_sig, 4)

writeData(wb, "S4a_N1_DE_genes",
          "Supplementary Table S4a. Differentially expressed genes (FDR < 0.05): N1 vs N0, adjusted for T stage, Gleason, and tissue source site",
          startRow = 1, startCol = 1)
addStyle(wb, "S4a_N1_DE_genes", title_style, rows = 1, cols = 1)
writeData(wb, "S4a_N1_DE_genes", de_n1_sig, startRow = 3, startCol = 1,
          headerStyle = hdr_style)

fg_n1 <- safe_read("tables/08b/FGSEA_Hallmark_N1_vs_N0_adjTGleason_TSS.csv")
fg_n1_sig <- fg_n1[fg_n1$padj < 0.05, ]
fg_n1_sig <- format_table(fg_n1_sig, 4)

writeData(wb, "S4b_N1_fgsea",
          "Supplementary Table S4b. Hallmark pathways enriched at FDR < 0.05: N1 vs N0 (adjusted)",
          startRow = 1, startCol = 1)
addStyle(wb, "S4b_N1_fgsea", title_style, rows = 1, cols = 1)
writeData(wb, "S4b_N1_fgsea", fg_n1_sig, startRow = 3, startCol = 1,
          headerStyle = hdr_style)

concord <- safe_read("tables/08b/Pathway_overlap_N1_vs_SPAG1.csv")
concord <- format_table(concord, 4)

writeData(wb, "S4c_Concordance",
          "Supplementary Table S4c. Pathway-level concordance between adjusted N1 vs N0 and SPAG1 Q4 vs Q1 analyses",
          startRow = 1, startCol = 1)
addStyle(wb, "S4c_Concordance", title_style, rows = 1, cols = 1)
writeData(wb, "S4c_Concordance", concord, startRow = 3, startCol = 1,
          headerStyle = hdr_style)



# =====================================================================
# Supp Table — Internal validation 
# =====================================================================
addWorksheet(wb, "S5_Internal_validation")

auc_models <- safe_read("tables/15/AUC_models.csv")
auc_delta  <- safe_read("tables/15/AUC_delta_and_pvalue.csv")
auc_boot   <- safe_read("tables/15/AUC_delta_bootstrap_CI.csv")
val_full   <- safe_read("tables/12_validation/Table_internal_validation.csv")

# Pull values from the wide-format validation table
auc_full_apparent  <- val_full$Apparent[val_full$Metric == "AUC (clinical + SPAG1)"]
auc_full_corrected <- val_full$Optimism_corrected[val_full$Metric == "AUC (clinical + SPAG1)"]
auc_full_optimism  <- val_full$Mean_optimism[val_full$Metric == "AUC (clinical + SPAG1)"]
auc_baseline       <- val_full$Apparent[val_full$Metric == "AUC (clinical only — m0)"]
calib_slope        <- val_full$Optimism_corrected[val_full$Metric == "Calibration slope"]

# Extract Hosmer-Lemeshow p from the Note column
hl_note <- val_full$Note[val_full$Metric == "Hosmer-Lemeshow test"]
hl_p    <- as.numeric(sub(".*p\\s*=\\s*([0-9.]+).*", "\\1", hl_note))

val_summary <- data.frame(
  Metric = c(
    "AUC, clinical baseline (T + Gleason)",
    "AUC, clinical + SPAG1 (full model)",
    "Delta AUC (DeLong test)",
    "Bootstrap 95% CI for Delta AUC",
    "DeLong p-value",
    "Optimism-corrected AUC (B = 2000)",
    "Optimism (apparent - corrected)",
    "Calibration slope (optimism-corrected)",
    "Hosmer-Lemeshow p-value"
  ),
  Value = c(
    sprintf("%.3f", auc_baseline),
    sprintf("%.3f", auc_full_apparent),
    sprintf("%.3f", auc_delta$Delta_AUC[1]),
    sprintf("%.3f - %.3f", auc_boot$CI_low[1], auc_boot$CI_high[1]),
    format(auc_delta$DeLong_p_value[1], digits = 3, scientific = TRUE),
    sprintf("%.3f", auc_full_corrected),
    sprintf("%.3f", auc_full_optimism),
    sprintf("%.3f", calib_slope),
    sprintf("%.3f", hl_p)
  ),
  stringsAsFactors = FALSE
)

writeData(wb, "S5_Internal_validation",
          "Supplementary Table. Discrimination and internal validation of the primary logistic regression model",
          startRow = 1, startCol = 1)
addStyle(wb, "S5_Internal_validation", title_style, rows = 1, cols = 1)
writeData(wb, "S5_Internal_validation", val_summary, startRow = 3, startCol = 1,
          headerStyle = hdr_style)


# =====================================================================
# Supp Table — Missingness
# =====================================================================
addWorksheet(wb, "S6_Missingness")
miss <- safe_read("tables/11a/Missingness_SPAG1_included_vs_excluded.csv")
miss <- format_table(miss, 3)

writeData(wb, "S6_Missingness",
          "Supplementary Table. Comparison of patients included in primary analysis (n = 420) versus those excluded for indeterminate nodal status (NX, n = 73)",
          startRow = 1, startCol = 1)
addStyle(wb, "S6_Missingness", title_style, rows = 1, cols = 1)
writeData(wb, "S6_Missingness", miss, startRow = 3, startCol = 1,
          headerStyle = hdr_style)

# =====================================================================
# Supp Table — Subgroup consistency
# =====================================================================
addWorksheet(wb, "S7_Subgroup_consistency")
sub <- safe_read("tables/11/Subgroup_consistency_SPAG1.csv")
sub <- format_table(sub, 3)

writeData(wb, "S7_Subgroup_consistency",
          "Supplementary Table. SPAG1 odds ratio for N1 within clinical subgroups",
          startRow = 1, startCol = 1)
addStyle(wb, "S7_Subgroup_consistency", title_style, rows = 1, cols = 1)
writeData(wb, "S7_Subgroup_consistency", sub, startRow = 3, startCol = 1,
          headerStyle = hdr_style)

# =====================================================================
# Supp Table S8 — Microenvironment correlations
# =====================================================================
addWorksheet(wb, "S8_Microenvironment")

# ESTIMATE correlations from Script 07
# TIMER2.0 purity correlation obtained from web portal
# https://timer.cistrome.org/ — accessed April 2026
micro_tbl <- data.frame(
  Variable     = c("ImmuneScore", "StromalScore", "Tumour purity"),
  Source       = c("ESTIMATE", "ESTIMATE", "TIMER2.0"),
  n            = c(497L, 497L, 497L),
  Spearman_rho = c(0.056, 0.030, -0.017),
  p_value      = c(0.212, 0.473, 0.711),
  stringsAsFactors = FALSE
)

writeData(wb, "S8_Microenvironment",
          "Supplementary Table S8. Spearman correlations between SPAG1 expression and tumour microenvironmental composition in TCGA-PRAD (n = 497). ImmuneScore and StromalScore were derived from the ESTIMATE algorithm. Tumour purity estimates were obtained from the TIMER2.0 web portal (https://timer.cistrome.org/).",
          startRow = 1, startCol = 1)
addStyle(wb, "S8_Microenvironment", title_style, rows = 1, cols = 1)
writeData(wb, "S8_Microenvironment", micro_tbl, startRow = 3, startCol = 1,
          headerStyle = hdr_style)
# =====================================================================
# Auto-size columns and save
# =====================================================================
for (sheet in names(wb)) {
  setColWidths(wb, sheet, cols = 1:20, widths = "auto")
}

output_path <- "tables/Supplementary_Tables.xlsx"
saveWorkbook(wb, output_path, overwrite = TRUE)
cat("\nSaved consolidated supplementary tables to:", output_path, "\n")
cat("Sheets created:", length(names(wb)), "\n")
cat("Sheet names:\n")
for (s in names(wb)) cat("  -", s, "\n")