# ============================================
# 11_SPAG1_logistic_regression_N1.R
# ============================================
# Purpose:
#   A) Primary logistic regression: SPAG1 predicts N1 beyond
#      clinical variables (T stage + Gleason grade)
#   B) Subgroup consistency: SPAG1-N1 association across
#      T stage and Gleason grade strata
#   C) Continuous vs categorical robustness: confirms
#      dose-response pattern and linear assumption validity
#   D) Model assumption checks
#
# Inputs:
#   data_processed/prad_master_features.rds
#
# Outputs:
#   tables/11/Logistic_N1_primary_model.csv
#   tables/11/N_stage_unadjusted.csv
#   tables/11/Logistic_N1_OR_table.csv
#   tables/11/Logistic_N1_model_coefficients_broom.csv
#   tables/11/Subgroup_consistency_SPAG1.csv
#   tables/11/SPAG1_quartile_robustness.csv
#   tables/11/SPAG1_quartile_N1_rates.csv
# ============================================

library(dplyr)
library(broom)
library(car)
library(splines)
library(detectseparation)

dir.create("tables/11", recursive = TRUE, showWarnings = FALSE)


# ----------------------------
# Load master dataset
# ----------------------------

dat <- readRDS("data_processed/prad_master_features.rds")

# ----------------------------
# A) Logistic regression for N1
# ----------------------------

# Make/refresh outcome
dat <- dat %>%
  mutate(
    N_binary = dplyr::case_when(
      tolower(as.character(N_group2)) == "n1" ~ 1,
      tolower(as.character(N_group2)) == "n0" ~ 0,
      TRUE ~ NA_real_
    )
  )

# ---- Unadjusted model: SPAG1 only ----
dat_N_unadj <- dat %>% filter(!is.na(N_binary), !is.na(SPAG1))
dat_N_unadj$N_group2 <- factor(dat_N_unadj$N_group2, levels = c("N0", "N1"))

m_N_unadj <- glm(N_group2 ~ SPAG1, data = dat_N_unadj, family = binomial())
res_N_unadj <- broom::tidy(m_N_unadj, exponentiate = TRUE, conf.int = TRUE)
write.csv(res_N_unadj, "tables/11/N_stage_unadjusted.csv", row.names = FALSE)
message("Saved: tables/11/N_stage_unadjusted.csv")


# Predictors
vars <- c("N_binary", "SPAG1", "T_group2", "Gleason_group2")
missing_vars <- setdiff(vars, names(dat))
if (length(missing_vars) > 0) stop("Missing columns in dat: ", paste(missing_vars, collapse = ", "))

# Complete-case dataset 
dat_n1 <- dat[complete.cases(dat[, vars]), , drop = FALSE]

cat("\n=== Logistic model dataset ===\n")
cat("Total samples in dat:", nrow(dat), "\n")
cat("Samples used in dat_n1 (complete cases):", nrow(dat_n1), "\n")
cat("N1 count:", sum(dat_n1$N_binary == 1), "\n")
cat("N0 count:", sum(dat_n1$N_binary == 0), "\n\n")

# ---- Primary model: clinical + SPAG1 ----
m0 <- glm(N_binary ~ T_group2 + Gleason_group2,
          data = dat_n1, family = binomial)

m1 <- glm(N_binary ~ T_group2 + Gleason_group2 + SPAG1,
          data = dat_n1, family = binomial)

fit_n1 <- m1

cat("\nPrimary model: T_group2 + Gleason_group2 + SPAG1\n")
print(summary(fit_n1))

# LRT: SPAG1 vs clinical baseline
lrt_m1 <- anova(m0, m1, test = "LRT")
cat("\nLRT — SPAG1 vs clinical baseline:\n")
print(lrt_m1)

# ---- Save full model coefficients (broom — profile likelihood CIs) ----
coeff_tbl <- broom::tidy(m1, exponentiate = TRUE, conf.int = TRUE) %>%
  select(term, estimate, conf.low, conf.high, p.value) %>%
  mutate(
    OR      = round(estimate, 3),
    CI_low  = round(conf.low, 3),
    CI_high = round(conf.high, 3),
    p_value = signif(p.value, 3)
  ) %>%
  select(term, OR, CI_low, CI_high, p_value)

write.csv(coeff_tbl,
          "tables/11/Logistic_N1_model_coefficients_broom.csv",
          row.names = FALSE)
message("Saved: tables/11/Logistic_N1_model_coefficients_broom.csv")

cat("\nFull model coefficients (profile likelihood CIs):\n")
print(coeff_tbl, row.names = FALSE)

# Save primary model summary
primary_tbl <- data.frame(
  model       = c("m0_clinical", "m1_clinical+SPAG1"),
  df          = c(m0$df.residual, m1$df.residual),
  AIC         = round(c(AIC(m0), AIC(m1)), 3),
  Delta_AIC   = round(c(0, AIC(m1) - AIC(m0)), 3),
  LRT_p_vs_m0 = c(NA, lrt_m1$`Pr(>Chi)`[2]),
  stringsAsFactors = FALSE
)

write.csv(primary_tbl,
          "tables/11/Logistic_N1_primary_model.csv",
          row.names = FALSE)
message("Saved: tables/11/Logistic_N1_primary_model.csv")
cat("\nModel comparison:\n")
print(primary_tbl, row.names = FALSE)

# =============================================================
# B) Subgroup consistency analysis
# =============================================================
# Purpose: assess whether SPAG1-N1 association is consistent
# across T stage and Gleason grade strata.
#
# Design:
#   - Same complete-case dataset as primary analysis (n=420)
#   - Unadjusted SPAG1-only model within each stratum
#   - T2 stratum excluded if N1 events < 10 (EPV < 10)
#
# Output:
#   tables/11/Subgroup_consistency_SPAG1.csv
# =============================================================

# Complete-case dataset — same as primary analysis
subgroup_vars <- c("N_binary", "T_group2", "Gleason_group2", "SPAG1")
dat_sub <- dat[complete.cases(dat[, subgroup_vars]), , drop = FALSE]

# Recreate N_binary as integer if needed for consistency
dat_sub <- dat_sub %>%
  mutate(
    N_binary = case_when(
      tolower(as.character(N_group2)) == "n1" ~ 1L,
      tolower(as.character(N_group2)) == "n0" ~ 0L,
      TRUE ~ NA_integer_
    ),
    T_group2       = factor(T_group2,
                            levels = c("T2", "T3_4")),
    Gleason_group2 = factor(Gleason_group2,
                            levels = c("Low_intermediate", "High"))
  )

cat("=============================================================\n")
cat(" SUBGROUP CONSISTENCY — SPAG1 vs N1\n")
cat("=============================================================\n")
cat(sprintf(" Full dataset: n=%d | N1=%d | N0=%d\n\n",
            nrow(dat_sub),
            sum(dat_sub$N_binary == 1),
            sum(dat_sub$N_binary == 0)))

# Helper function
run_subgroup <- function(data, label) {
  n    <- nrow(data)
  n_N1 <- sum(data$N_binary == 1)
  n_N0 <- sum(data$N_binary == 0)
  epv  <- n_N1 / 1
  
  cat(sprintf(" --- %s ---\n", label))
  cat(sprintf("  n=%d | N1=%d | N0=%d | EPV=%.1f\n",
              n, n_N1, n_N0, epv))
  
  if (n_N1 < 10) {
    cat(sprintf(
      "  SKIP: only %d N1 events (EPV=%.1f) — insufficient for stable estimates\n\n",
      n_N1, epv))
    return(data.frame(
      Stratum   = label,
      n         = n,
      N1        = n_N1,
      EPV       = round(epv, 1),
      OR        = NA_real_,
      CI_low    = NA_real_,
      CI_high   = NA_real_,
      p         = NA_real_,
      Note      = "Insufficient N1 events",
      stringsAsFactors = FALSE
    ))
  }
  
  fit <- tryCatch(
    glm(N_binary ~ SPAG1, data = data, family = binomial()),
    error = function(e) {
      cat("  ERROR:", conditionMessage(e), "\n\n")
      return(NULL)
    }
  )
  
  if (is.null(fit)) return(NULL)
  
  res <- broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == "SPAG1") %>%
    transmute(
      Stratum = label,
      n       = n,
      N1      = n_N1,
      EPV     = round(epv, 1),
      OR      = round(estimate, 3),
      CI_low  = round(conf.low, 3),
      CI_high = round(conf.high, 3),
      p       = signif(p.value, 3),
      Note    = ""
    )
  
  cat(sprintf("  OR = %.3f (95%% CI %.3f-%.3f)  p = %s\n\n",
              res$OR, res$CI_low, res$CI_high,
              format(res$p, scientific = TRUE)))
  return(as.data.frame(res))
}

# Run four strata
subgroup_results <- bind_rows(
  run_subgroup(dat_sub %>% filter(T_group2 == "T2"),
               "T2 only"),
  run_subgroup(dat_sub %>% filter(T_group2 == "T3_4"),
               "T3/T4 only"),
  run_subgroup(dat_sub %>% filter(Gleason_group2 == "Low_intermediate"),
               "Low/Intermediate Gleason"),
  run_subgroup(dat_sub %>% filter(Gleason_group2 == "High"),
               "High Gleason")
)

# Add primary model as reference row
primary_ref <- broom::tidy(m1, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term == "SPAG1") %>%
  transmute(
    Stratum = "Full cohort (adjusted: T + Gleason + SPAG1)",
    n       = nrow(dat_sub),
    N1      = sum(dat_sub$N_binary == 1),
    EPV     = round(sum(dat_sub$N_binary == 1) / 3, 1),
    OR      = round(estimate, 3),
    CI_low  = round(conf.low, 3),
    CI_high = round(conf.high, 3),
    p       = signif(p.value, 3),
    Note    = "Primary adjusted model — reference"
  )

final_tbl <- bind_rows(as.data.frame(primary_ref), subgroup_results)

write.csv(final_tbl,
          "tables/11/Subgroup_consistency_SPAG1.csv",
          row.names = FALSE)
message("Saved: tables/11/Subgroup_consistency_SPAG1.csv")

cat("=============================================================\n")
cat(" SUBGROUP SUMMARY\n")
cat("=============================================================\n")
print(final_tbl[, c("Stratum", "n", "N1", "EPV",
                    "OR", "CI_low", "CI_high", "p", "Note")],
      row.names = FALSE)


# =============================================================
# Section C: Continuous vs categorical robustness 
# =============================================================
# Purpose: confirm SPAG1-N1 association follows a dose-response
# pattern and that the continuous model assumption is valid.
#
# Three complementary specifications:
#   1. Continuous SPAG1 — primary model
#   2. Quartile categories — Q1 as reference, tests non-linearity
#   3. Trend test — ordered integer scores, tests dose-response
#
# Quartiles defined within the 420-patient analytic sample
# to avoid information leakage from full dataset.
#
# Output:
#   tables/11/SPAG1_quartile_robustness.csv
# =============================================================

# Define quartiles within analytic sample
dat_n1 <- dat_n1 %>%
  mutate(
    SPAG1_Q        = ntile(SPAG1, 4),
    SPAG1_Q_factor = factor(SPAG1_Q,
                            levels = 1:4,
                            labels = c("Q1", "Q2", "Q3", "Q4"))
  )

cat("=============================================================\n")
cat(" CONTINUOUS VS QUARTILE ROBUSTNESS\n")
cat("=============================================================\n")

# Raw N1 rate per quartile
quartile_summary <- dat_n1 %>%
  group_by(SPAG1_Q_factor) %>%
  summarise(
    n      = n(),
    N1     = sum(N_binary == 1),
    N0     = sum(N_binary == 0),
    N1_pct = round(100 * mean(N_binary == 1), 1),
    .groups = "drop"
  )

cat("\n SPAG1 quartile boundaries (analytic sample n=420):\n")
print(quantile(dat_n1$SPAG1, probs = seq(0, 1, 0.25)))
cat("\n N1 rate per quartile:\n")
print(quartile_summary, row.names = FALSE)

# Model 1: continuous (reference — same as m2)
m_cont <- m1

# Model 2: quartile categories
m_cat <- glm(N_binary ~ T_group2 + Gleason_group2 + SPAG1_Q_factor,
             data = dat_n1, family = binomial())

# Model 3: trend test
m_trend <- glm(N_binary ~ T_group2 + Gleason_group2 + SPAG1_Q,
               data = dat_n1, family = binomial())

# Extract results
cont_res <- tidy(m_cont, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term == "SPAG1") %>%
  transmute(
    Specification = "Continuous (per 1 log2-RSEM unit)",
    OR      = round(estimate, 3),
    CI_low  = round(conf.low, 3),
    CI_high = round(conf.high, 3),
    p       = signif(p.value, 3),
    AIC     = round(AIC(m_cont), 2)
  )

cat_res <- tidy(m_cat, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(grepl("SPAG1", term)) %>%
  transmute(
    Specification = case_match(term,
                               "SPAG1_Q_factorQ2" ~ "Q2 vs Q1",
                               "SPAG1_Q_factorQ3" ~ "Q3 vs Q1",
                               "SPAG1_Q_factorQ4" ~ "Q4 vs Q1"
    ),
    OR      = round(estimate, 3),
    CI_low  = round(conf.low, 3),
    CI_high = round(conf.high, 3),
    p       = signif(p.value, 3),
    AIC     = round(AIC(m_cat), 2)
  )

trend_res <- tidy(m_trend, exponentiate = TRUE, conf.int = TRUE) %>%
  filter(term == "SPAG1_Q") %>%
  transmute(
    Specification = "Trend (per quartile increment)",
    OR      = round(estimate, 3),
    CI_low  = round(conf.low, 3),
    CI_high = round(conf.high, 3),
    p       = signif(p.value, 3),
    AIC     = round(AIC(m_trend), 2)
  )

robustness_tbl <- bind_rows(cont_res, cat_res, trend_res)

cat("\n Model results:\n")
print(robustness_tbl, row.names = FALSE)

# AIC comparison
cat(sprintf("\n AIC — Continuous: %.2f | Quartile: %.2f | Trend: %.2f\n",
            AIC(m_cont), AIC(m_cat), AIC(m_trend)))

# Interpretation
cat("\n Interpretation:\n")
cat(sprintf(" Q1 N1 rate: %.1f%%  →  Q4 N1 rate: %.1f%%\n",
            quartile_summary$N1_pct[1],
            quartile_summary$N1_pct[4]))
cat(sprintf(" Trend OR per quartile increment: %.3f (p=%s)\n",
            trend_res$OR,
            format(trend_res$p, scientific = TRUE)))

# Save
write.csv(robustness_tbl,
          "tables/11/SPAG1_quartile_robustness.csv",
          row.names = FALSE)
write.csv(quartile_summary,
          "tables/11/SPAG1_quartile_N1_rates.csv",
          row.names = FALSE)
message("Saved: tables/11/SPAG1_quartile_robustness.csv")
message("Saved: tables/11/SPAG1_quartile_N1_rates.csv")

# ---- Combined OR table: unadjusted and adjusted SPAG1 ----
# Uses broom profile likelihood CIs — more accurate than Wald
# especially for sparse cells (T2 with only 3 N1 events)
extract_spag1_or <- function(model, label) {
  broom::tidy(model, exponentiate = TRUE, conf.int = TRUE) %>%
    filter(term == "SPAG1") %>%
    transmute(
      Model   = label,
      OR      = round(estimate, 3),
      CI_low  = round(conf.low, 3),
      CI_high = round(conf.high, 3),
      p_value = signif(p.value, 3)
    ) %>%
    as.data.frame()
}


or_combined <- bind_rows(
  extract_spag1_or(m_N_unadj,
                   "Unadjusted (SPAG1 only)"),
  extract_spag1_or(m1,
                   "Adjusted (T stage + Gleason + SPAG1)")
)

write.csv(or_combined,
          "tables/11/SPAG1_OR_all_models.csv",
          row.names = FALSE)
message("Saved: tables/11/SPAG1_OR_all_models.csv")
cat("\nSPAG1 OR across model specifications:\n")
print(or_combined, row.names = FALSE)

# =============================================================
# D) Model assumption checks
# =============================================================
# Purpose: confirm logistic regression assumptions are met.
# Results reported in Methods as a single summary statement.
# Checks: multicollinearity (VIF), convergence, separation,
#         linearity of SPAG1 (spline vs linear LRT + AIC)
# =============================================================

cat("=============================================================\n")
cat(" MODEL ASSUMPTION CHECKS\n")
cat("=============================================================\n")

# ---- 1. Convergence ----
cat(sprintf("\n1. Convergence: %s\n", 
            ifelse(m1$converged, "TRUE — no issues", "FALSE — INVESTIGATE")))

# ---- 2. Multicollinearity (VIF) ----
vif_vals <- car::vif(m1)
cat("\n2. Multicollinearity (VIF):\n")
print(round(vif_vals, 3))
if (any(vif_vals > 5)) {
  warning("VIF > 5 detected — collinearity concern")
} else {
  cat(sprintf("   Max VIF = %.3f — no collinearity concern\n", max(vif_vals)))
}

# ---- 3. Separation check ----
sep_check <- glm(
  N_binary ~ T_group2 + Gleason_group2 + SPAG1,
  data   = dat_n1,
  family = binomial(),
  method = detectseparation::detect_separation
)
cat(sprintf("\n3. Separation: %s\n",
            ifelse(sep_check$separation, 
                   "TRUE — SEPARATION DETECTED", 
                   "FALSE — no separation")))

# ---- 4. Linearity of SPAG1 (spline vs linear LRT) ----
m_spline <- glm(
  N_binary ~ T_group2 + Gleason_group2 + splines::ns(SPAG1, df = 3),
  data   = dat_n1,
  family = binomial()
)
spline_lrt <- anova(m1, m_spline, test = "Chisq")
spline_p   <- spline_lrt$`Pr(>Chi)`[2]
spline_dAIC <- AIC(m_spline) - AIC(m1)

cat(sprintf(
  "\n4. Linearity of SPAG1 — spline vs linear:\n   LRT p = %.4f | Delta AIC = %.2f\n",
  spline_p, spline_dAIC
))
if (spline_p > 0.05 && spline_dAIC > -2) {
  cat("   Linear term adequate — spline does not improve fit\n")
} else {
  cat("   WARNING: spline may improve fit — consider non-linear term\n")
}

# ---- 5. Predicted probabilities ----
pred_prob <- fitted(m1)
cat(sprintf(
  "\n5. Predicted probabilities: min=%.4f, max=%.4f\n",
  min(pred_prob), max(pred_prob)
))
if (sum(pred_prob < 0.001) > 0 || sum(pred_prob > 0.999) > 0) {
  cat("   WARNING: extreme predicted probabilities detected\n")
} else {
  cat("   No boundary predictions — no concern\n")
}

cat("\n=============================================================\n")
cat(" Assumption checks complete. No major violations expected.\n")
cat("=============================================================\n")
cat("\n=============================================================\n")
cat(" Script 11 complete.\n")
cat("=============================================================\n")