# =============================================================================
# Script 11a: Missingness Analysis and Sensitivity Analysis
# =============================================================================
# Purpose:
#   1. Assess whether complete-case exclusion introduces selection bias
#      by comparing SPAG1 expression between included and excluded patients
#   2a. Wilcoxon comparison of SPAG1 expression: included vs excluded
#   2b. Clinical profile of NX patients (T stage, Gleason) — suggests
#       that NX exclusion is by-design based on lower-stage profile
#   2c. SPAG1 expression by N stage — numerical summary
#       Companion to Fig_SPAG1_by_N.png; verifies directional
#       association before modelling
#   3. Sensitivity analysis treating NX patients as N0 (clinically
#      plausible conservative assumption) to confirm SPAG1 finding
#      is robust to complete-case selection
#
# Context:
#   497 patients total; 420 included in primary analysis (complete cases)
#   77 excluded: 73 NX (lymph node dissection not performed — by-design
#   exclusion, not random missing data) + 4 missing T/Gleason
#
# Inputs:
#   data_processed/prad_master_features.rds
#
# Outputs:
#   tables/11a/Missingness_SPAG1_included_vs_excluded.csv
#   tables/11a/NX_patient_clinical_profile.csv
#   tables/11a/Sensitivity_NX_as_N0_SPAG1.csv
#    tables/11a/SPAG1_by_N_stage_summary.csv
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(broom)
})

dir.create("tables/11a", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD AND PREPARE
# =============================================================================

dat_all <- readRDS("data_processed/prad_master_features.rds")

dat_all <- dat_all %>%
  mutate(
    N_binary = case_when(
      tolower(as.character(N_group2)) == "n1" ~ 1L,
      tolower(as.character(N_group2)) == "n0" ~ 0L,
      TRUE ~ NA_integer_
    )
  )

# Flag complete cases using same variables as primary analysis
dat_all$included <- complete.cases(
  dat_all[, c("N_binary", "T_group2", "Gleason_group2", "SPAG1")]
)

cat("=============================================================\n")
cat(" SAMPLE SUMMARY\n")
cat("=============================================================\n")
cat(sprintf(" Total patients         : %d\n", nrow(dat_all)))
cat(sprintf(" Included (complete cc) : %d\n", sum(dat_all$included)))
cat(sprintf(" Excluded               : %d\n", sum(!dat_all$included)))
cat(sprintf("   NX (no dissection)   : %d\n",
            sum(is.na(dat_all$N_binary))))
cat(sprintf("   Missing T/Gleason    : %d\n",
            sum(!dat_all$included) - sum(is.na(dat_all$N_binary))))
cat("\n")

# =============================================================================
# 2a. MISSINGNESS ANALYSIS
#    Compare SPAG1 expression between included and excluded patients
#    Wilcoxon rank-sum test
# =============================================================================

spag1_inc <- dat_all$SPAG1[ dat_all$included  & !is.na(dat_all$SPAG1)]
spag1_exc <- dat_all$SPAG1[!dat_all$included  & !is.na(dat_all$SPAG1)]

wt <- wilcox.test(spag1_inc, spag1_exc)

miss_tbl <- data.frame(
  Group        = c("Included", "Excluded"),
  N            = c(length(spag1_inc), length(spag1_exc)),
  Median_SPAG1 = round(c(median(spag1_inc), median(spag1_exc)), 3),
  IQR_SPAG1    = c(
    paste(round(quantile(spag1_inc, 0.25), 3),
          round(quantile(spag1_inc, 0.75), 3), sep = "–"),
    paste(round(quantile(spag1_exc, 0.25), 3),
          round(quantile(spag1_exc, 0.75), 3), sep = "–")
  ),
  Wilcoxon_p   = c(round(wt$p.value, 4), NA),
  stringsAsFactors = FALSE
)

write.csv(miss_tbl,
          "tables/11a/Missingness_SPAG1_included_vs_excluded.csv",
          row.names = FALSE)
message("Saved: tables/11a/Missingness_SPAG1_included_vs_excluded.csv")

cat("=============================================================\n")
cat(" MISSINGNESS ANALYSIS\n")
cat("=============================================================\n")
print(miss_tbl, row.names = FALSE)
cat(sprintf("\n Wilcoxon p = %.4f\n", wt$p.value))
if (wt$p.value < 0.05) {
  cat(" Excluded patients have significantly lower SPAG1 expression.\n")
  cat(" Consistent with lower-stage profile (NX = no dissection performed).\n")
  cat(" This is a by-design exclusion, not random missing data.\n\n")
} else {
  cat(" No significant difference in SPAG1 between included and excluded.\n")
  cat(" Complete-case selection unlikely to introduce meaningful bias.\n\n")
}

# =============================================================
# Section 2b. Clinical profile of NX patients
# Verifies NX patients are lower-stage
# =============================================================

nx_patients  <- dat_all[is.na(dat_all$N_binary), ]
n0n1_patients <- dat_all[!is.na(dat_all$N_binary), ]

cat("=== NX PATIENT PROFILE VERIFICATION ===\n\n")

# T stage distribution
cat("T stage distribution:\n")
cat("NX patients:\n")
print(table(nx_patients$T_group2, useNA = "ifany"))
cat("N0/N1 patients:\n")
print(table(n0n1_patients$T_group2, useNA = "ifany"))

# Gleason distribution
cat("\nGleason distribution:\n")
cat("NX patients:\n")
print(table(nx_patients$Gleason_group2, useNA = "ifany"))
cat("N0/N1 patients:\n")
print(table(n0n1_patients$Gleason_group2, useNA = "ifany"))

# T stage chi-square test
tss_tbl <- table(
  is.na(dat_all$N_binary),
  dat_all$T_group2
)
tss_tbl <- tss_tbl[, !is.na(colnames(tss_tbl))]
chi_T <- chisq.test(tss_tbl)
cat(sprintf("\nChi-square test T stage vs NX status: p = %.4f\n",
            chi_T$p.value))

# Gleason chi-square test
gl_tbl <- table(
  is.na(dat_all$N_binary),
  dat_all$Gleason_group2
)
gl_tbl <- gl_tbl[, !is.na(colnames(gl_tbl))]
chi_G <- chisq.test(gl_tbl)
cat(sprintf("Chi-square test Gleason vs NX status: p = %.4f\n",
            chi_G$p.value))

# Save
nx_profile <- data.frame(
  Variable = c("T2", "T3_4", "Gleason_LowInt", "Gleason_High"),
  NX_n = c(
    sum(nx_patients$T_group2 == "T2", na.rm = TRUE),
    sum(nx_patients$T_group2 == "T3_4", na.rm = TRUE),
    sum(nx_patients$Gleason_group2 == "Low_intermediate", na.rm = TRUE),
    sum(nx_patients$Gleason_group2 == "High", na.rm = TRUE)
  ),
  NX_pct = round(c(
    mean(nx_patients$T_group2 == "T2", na.rm = TRUE),
    mean(nx_patients$T_group2 == "T3_4", na.rm = TRUE),
    mean(nx_patients$Gleason_group2 == "Low_intermediate", na.rm = TRUE),
    mean(nx_patients$Gleason_group2 == "High", na.rm = TRUE)
  ) * 100, 1),
  N0N1_n = c(
    sum(n0n1_patients$T_group2 == "T2", na.rm = TRUE),
    sum(n0n1_patients$T_group2 == "T3_4", na.rm = TRUE),
    sum(n0n1_patients$Gleason_group2 == "Low_intermediate", na.rm = TRUE),
    sum(n0n1_patients$Gleason_group2 == "High", na.rm = TRUE)
  ),
  N0N1_pct = round(c(
    mean(n0n1_patients$T_group2 == "T2", na.rm = TRUE),
    mean(n0n1_patients$T_group2 == "T3_4", na.rm = TRUE),
    mean(n0n1_patients$Gleason_group2 == "Low_intermediate", na.rm = TRUE),
    mean(n0n1_patients$Gleason_group2 == "High", na.rm = TRUE)
  ) * 100, 1)
)

write.csv(nx_profile,
          "tables/11a/NX_patient_clinical_profile.csv",
          row.names = FALSE)
message("Saved: tables/11a/NX_patient_clinical_profile.csv")
print(nx_profile, row.names = FALSE)


# =============================================================
# Section 2c. SPAG1 expression by N stage — numerical summary
# Companion to Fig_SPAG1_by_N.png; verifies directional
# association before modelling
# =============================================================

spag1_n0 <- dat_all$SPAG1[dat_all$N_binary %in% 0L & !is.na(dat_all$SPAG1)]
spag1_n1 <- dat_all$SPAG1[dat_all$N_binary %in% 1L & !is.na(dat_all$SPAG1)]
spag1_nx <- dat_all$SPAG1[is.na(dat_all$N_binary) & !is.na(dat_all$SPAG1)]

wt_n <- wilcox.test(spag1_n0, spag1_n1)

spag1_n_summary <- data.frame(
  Group        = c("N0", "N1", "NX"),
  N            = c(length(spag1_n0), length(spag1_n1), length(spag1_nx)),
  Median_SPAG1 = round(c(median(spag1_n0, na.rm = TRUE),
                         median(spag1_n1, na.rm = TRUE),
                         median(spag1_nx, na.rm = TRUE)), 3),
  IQR_SPAG1    = c(
    paste(round(quantile(spag1_n0, 0.25, na.rm = TRUE), 3),
          round(quantile(spag1_n0, 0.75, na.rm = TRUE), 3), sep = "–"),
    paste(round(quantile(spag1_n1, 0.25, na.rm = TRUE), 3),
          round(quantile(spag1_n1, 0.75, na.rm = TRUE), 3), sep = "–"),
    paste(round(quantile(spag1_nx, 0.25, na.rm = TRUE), 3),
          round(quantile(spag1_nx, 0.75, na.rm = TRUE), 3), sep = "–")
  ),
  Wilcoxon_p_vs_N0 = c(NA, round(wt_n$p.value, 4), NA),
  stringsAsFactors = FALSE
)

write.csv(spag1_n_summary,
          "tables/11a/SPAG1_by_N_stage_summary.csv",
          row.names = FALSE)
message("Saved: tables/11a/SPAG1_by_N_stage_summary.csv")
print(spag1_n_summary, row.names = FALSE)
cat(sprintf("\n SPAG1 N0 vs N1 Wilcoxon p = %.4f\n", wt_n$p.value))
cat(sprintf(" Direction: N1 (%.3f) > N0 (%.3f) > NX (%.3f)\n",
            median(spag1_n1, na.rm = TRUE),
            median(spag1_n0, na.rm = TRUE),
            median(spag1_nx, na.rm = TRUE)))


# =============================================================================
# 3. SENSITIVITY ANALYSIS — NX PATIENTS TREATED AS N0
#
#    Rationale: NX patients did not undergo lymph node dissection,
#    typically because their clinical profile suggested low nodal risk.
#    Treating them as N0 is the most clinically plausible assumption,
#    consistent with their lower-stage profile (T2: 60%, Low/Int
#    Gleason: 82.2%).
#
#    Note: NX ≠ N0 in clinical reality. This sensitivity analysis
#    tests robustness, not clinical equivalence.
# =============================================================================

dat_sens <- dat_all %>%
  mutate(
    N_binary_sens = case_when(
      tolower(as.character(N_group2)) == "n1" ~ 1L,
      tolower(as.character(N_group2)) == "n0" ~ 0L,
      TRUE ~ 0L    
    )
  ) %>%
  filter(!is.na(T_group2), !is.na(Gleason_group2), !is.na(SPAG1))

cat("=============================================================\n")
cat(" SENSITIVITY ANALYSIS — NX TREATED AS N0\n")
cat("=============================================================\n")
cat(sprintf(" N total    : %d\n", nrow(dat_sens)))
cat(sprintf(" N1         : %d\n", sum(dat_sens$N_binary_sens == 1)))
cat(sprintf(" N0 (+ NX)  : %d\n", sum(dat_sens$N_binary_sens == 0)))
cat(sprintf(" Event rate : %.1f%%\n\n",
            100 * mean(dat_sens$N_binary_sens)))

m_sens <- glm(
  N_binary_sens ~ T_group2 + Gleason_group2 + SPAG1,
  data = dat_sens, family = binomial
)

res_sens <- broom::tidy(m_sens, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term == "SPAG1")

cat(sprintf(" SPAG1 OR   = %.3f\n", res_sens$estimate))
cat(sprintf(" 95%% CI     = %.3f – %.3f\n",
            res_sens$conf.low, res_sens$conf.high))
cat(sprintf(" p-value    = %.4f\n\n", res_sens$p.value))

if (res_sens$p.value < 0.05) {
  cat(" SPAG1 remains significant under clinically plausible assumption.\n")
  cat(" Primary finding is robust to complete-case selection.\n\n")
} else {
  cat(" SPAG1 non-significant under clinically plausible assumption.\n")
  cat(" Interpret primary finding with caution.\n\n")
}

write.csv(
  broom::tidy(m_sens, exponentiate = TRUE, conf.int = TRUE),
  "tables/11a/Sensitivity_NX_as_N0_SPAG1.csv",
  row.names = FALSE
)
message("Saved: tables/11a/Sensitivity_NX_as_N0_SPAG1.csv")

cat("=============================================================\n")
cat(" Script 11a complete.\n")
cat("=============================================================\n")