# =============================================================================
# Script 12: Internal Validation — Bootstrap Optimism Correction
# =============================================================================
# Study:      TCGA-PRAD RNA-seq cohort (Illumina HiSeqV2)
# Outcome:    Lymph node involvement (N1 = 1, N0 = 0)
# Model:      N_binary ~ T_group2 + Gleason_group2 + SPAG1
# Method:     Harrell's bootstrap optimism correction (B = 2000)
#
# Metrics corrected:
#   (1) AUC               — discrimination
#   (2) Calibration slope — overfitting / sharpness of predictions
#
# Baseline model (m0: T_group2 + Gleason_group2) bootstrapped in parallel
# to derive optimism-corrected ΔAUC (full model vs clinical baseline).
#
# Outputs include calibration plots (apparent and corrected) and a
# combined figure for manuscript submission.
#
# Calibration intercept — NOT bootstrapped.
#   Rationale: For any logistic regression fitted with an intercept,
#   the calibration intercept on the training data is identically 0 by
#   construction (mean predicted probability = observed event rate).
#   Bootstrapping this quantity produces numerical instability (near-
#   separation in some resamples → extreme offset estimates) without
#   adding interpretable information. The apparent value (0.00) is
#   reported directly as evidence of calibration-in-the-large.
#
# Goodness-of-fit:
#   Hosmer-Lemeshow test (g = 10 groups) is reported alongside the
#   calibration plot as required by many clinical journals.
#   Limitation: HL test has reduced power with small samples and is
#   sensitive to grouping choice; it is reported for completeness and
#   the calibration plot is the primary calibration assessment.
#
# Inputs:
#   data_processed/prad_master_features.rds
#
# Outputs:
#   tables/12_validation/Table_internal_validation.csv
#   tables/12_validation/HL_observed_vs_expected.csv
#   figures/12/Figure_calibration_plot_main_model.png
#   figures/12/Figure_calibration_plot_corrected.png
#   figures/12/Figure_5_calibration_combined.png
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(pROC)
  library(ggplot2)
  library(scales)
  library(ResourceSelection)   # hoslem.test()
  library(patchwork)
})

# =============================================================================
# 0. SETTINGS
# =============================================================================

SEED <- 123
B    <- 2000
HL_G <- 10        # Number of groups for Hosmer-Lemeshow test


dir.create("tables/12_validation", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/12",              recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. LOAD & PREPARE DATA
# =============================================================================

dat <- readRDS("data_processed/prad_master_features.rds")

dat <- dat %>%
  mutate(
    N_binary = case_when(
      tolower(as.character(N_group2)) == "n1" ~ 1L,
      tolower(as.character(N_group2)) == "n0" ~ 0L,
      TRUE                                    ~ NA_integer_
    )
  )

model_vars <- c("N_binary", "T_group2", "Gleason_group2", "SPAG1")


dat_model <- dat[complete.cases(dat[, model_vars]), , drop = FALSE]
n         <- nrow(dat_model)

cat("=============================================================\n")
cat(" SAMPLE SUMMARY\n")
cat("=============================================================\n")
cat(sprintf(" Total samples  : %d\n",   n))
cat(sprintf(" Events  (N1)   : %d\n",   sum(dat_model$N_binary == 1)))
cat(sprintf(" Controls (N0)  : %d\n",   sum(dat_model$N_binary == 0)))
cat(sprintf(" Event rate     : %.1f%%\n\n", mean(dat_model$N_binary) * 100))

# =============================================================================
# 2. FINAL MODEL
# =============================================================================

formula_main <- N_binary ~ T_group2 + Gleason_group2 + SPAG1

fit_final <- suppressWarnings(
  glm(formula_main, data = dat_model, family = binomial)
)

# ---- Baseline clinical model (m0) for comparison ----
formula_m0  <- N_binary ~ T_group2 + Gleason_group2
fit_m0      <- suppressWarnings(glm(formula_m0, data = dat_model, family = binomial))
prob_m0     <- predict(fit_m0, type = "response")
lp_m0       <- predict(fit_m0, type = "link")
auc_m0_apparent <- as.numeric(auc(roc(dat_model$N_binary, prob_m0, quiet = TRUE)))
cat(sprintf(" AUC clinical model (m0) apparent: %.4f\n", auc_m0_apparent))

lp_final   <- predict(fit_final, type = "link")
prob_final <- predict(fit_final, type = "response")

# =============================================================================
# 3. APPARENT PERFORMANCE
# =============================================================================

auc_apparent <- as.numeric(
  auc(roc(dat_model$N_binary, prob_final, quiet = TRUE))
)

cal_int_apparent <- suppressWarnings(
  coef(
    glm(N_binary ~ offset(lp_final), data = dat_model, family = binomial)
  )[["(Intercept)"]]
)

cal_slope_apparent <- suppressWarnings(
  coef(
    glm(N_binary ~ lp_final, data = dat_model, family = binomial)
  )[["lp_final"]]
)

cat("=============================================================\n")
cat(" APPARENT PERFORMANCE  (optimistic — before correction)\n")
cat("=============================================================\n")
cat(sprintf(" AUC                   : %.4f\n",   auc_apparent))
cat(sprintf(" Calibration intercept : %.4f  (ideal = 0)\n", cal_int_apparent))
cat(sprintf(" Calibration slope     : %.4f  (ideal = 1)\n\n", cal_slope_apparent))

# =============================================================================
# 4. HOSMER-LEMESHOW GOODNESS-OF-FIT TEST
#
#   H0: No significant difference between observed and predicted event rates
#       across HL_G quantile-based groups (i.e., model is well calibrated).
#   H1: Significant lack of fit.
#
#   Interpretation:
#     p > 0.05 → no evidence of poor calibration (desired result)
#     p < 0.05 → evidence of miscalibration
#
#   Important caveat: HL test has reduced power in small samples (<200 events)
#   and is sensitive to the choice of g. Result should be interpreted
#   alongside the calibration plot, not in isolation.
# =============================================================================

hl_test <- hoslem.test(dat_model$N_binary, prob_final, g = HL_G)

cat("=============================================================\n")
cat(sprintf(" HOSMER-LEMESHOW TEST  (g = %d groups)\n", HL_G))
cat("=============================================================\n")
cat(sprintf(" Chi-squared statistic : %.4f\n", hl_test$statistic))
cat(sprintf(" Degrees of freedom    : %d\n",   hl_test$parameter))
cat(sprintf(" p-value               : %.4f\n", hl_test$p.value))

if (hl_test$p.value >= 0.05) {
  cat(" Interpretation        : No significant lack of fit (p >= 0.05)\n\n")
} else {
  cat(" Interpretation        : Evidence of lack of fit (p < 0.05)\n\n")
}

# Observed vs expected table (supplementary material)
hl_table <- as.data.frame(cbind(hl_test$observed, hl_test$expected))
colnames(hl_table) <- c("Obs_N0", "Obs_N1", "Exp_N0", "Exp_N1")
hl_table$Group <- seq_len(nrow(hl_table))
hl_table <- hl_table[, c("Group", "Obs_N0", "Obs_N1", "Exp_N0", "Exp_N1")]
hl_table[, c("Exp_N0", "Exp_N1")] <- round(hl_table[, c("Exp_N0", "Exp_N1")], 2)

write.csv(
  hl_table,
  "tables/12_validation/HL_observed_vs_expected.csv",
  row.names = FALSE
)
cat(" HL observed vs expected saved: tables/12_validation/HL_observed_vs_expected.csv\n\n")

# =============================================================================
# 5. BOOTSTRAP OPTIMISM CORRECTION
# =============================================================================

opt_auc   <- numeric(B)
opt_slope <- numeric(B)
opt_auc_m0 <- numeric(B)
n_failed  <- 0L

cat("=============================================================\n")
cat(sprintf(" BOOTSTRAP OPTIMISM CORRECTION  (B = %d, seed = %d)\n", B, SEED))
cat("=============================================================\n")
set.seed(SEED)
for (b in seq_len(B)) {
  
  idx      <- sample.int(n, size = n, replace = TRUE)
  dat_boot <- dat_model[idx, , drop = FALSE]
  
  if (length(unique(dat_boot$N_binary)) < 2L) {
    n_failed     <- n_failed + 1L
    opt_auc[b]   <- NA_real_
    opt_slope[b] <- NA_real_
    next
  }
  
  fit_b <- tryCatch(
    suppressWarnings(glm(formula_main, data = dat_boot, family = binomial)),
    error = function(e) NULL
  )
  
  if (is.null(fit_b)) {
    n_failed     <- n_failed + 1L
    opt_auc[b]   <- NA_real_
    opt_slope[b] <- NA_real_
    next
  }
  
  lp_train <- predict(fit_b, newdata = dat_boot,  type = "link")
  lp_test  <- predict(fit_b, newdata = dat_model, type = "link")
  
  auc_train <- tryCatch(
    as.numeric(auc(roc(dat_boot$N_binary,  plogis(lp_train), quiet = TRUE))),
    error = function(e) NA_real_
  )
  auc_test <- tryCatch(
    as.numeric(auc(roc(dat_model$N_binary, plogis(lp_test),  quiet = TRUE))),
    error = function(e) NA_real_
  )
  
  slope_train <- tryCatch(
    suppressWarnings(
      coef(glm(N_binary ~ lp_train, data = dat_boot,  family = binomial))[2]
    ),
    error = function(e) NA_real_
  )
  slope_test <- tryCatch(
    suppressWarnings(
      coef(glm(N_binary ~ lp_test,  data = dat_model, family = binomial))[2]
    ),
    error = function(e) NA_real_
  )
  
  opt_auc[b]   <- auc_train   - auc_test
  opt_slope[b] <- slope_train - slope_test
  
  # Optimism for m0 (clinical baseline)
  fit_m0_b <- tryCatch(
    suppressWarnings(glm(formula_m0, data = dat_boot, family = binomial)),
    error = function(e) NULL
  )
  if (!is.null(fit_m0_b)) {
    auc_m0_train <- tryCatch(
      as.numeric(auc(roc(dat_boot$N_binary,  plogis(predict(fit_m0_b, newdata = dat_boot,  type = "link")), quiet = TRUE))),
      error = function(e) NA_real_
    )
    auc_m0_test <- tryCatch(
      as.numeric(auc(roc(dat_model$N_binary, plogis(predict(fit_m0_b, newdata = dat_model, type = "link")), quiet = TRUE))),
      error = function(e) NA_real_
    )
    opt_auc_m0[b] <- auc_m0_train - auc_m0_test
  }
  
  if (b %% 500 == 0) {
    cat(sprintf("  %4d / %d resamples completed\n", b, B))
  }
}

cat(sprintf(
  "\n  Done. Failed resamples: %d / %d (%.1f%%)\n\n",
  n_failed, B, 100 * n_failed / B
))


# =============================================================================
# 6. CORRECTED ESTIMATES
# =============================================================================

mean_opt_auc   <- mean(opt_auc,   na.rm = TRUE)
mean_opt_slope <- mean(opt_slope, na.rm = TRUE)
mean_opt_auc_m0 <- mean(opt_auc_m0, na.rm = TRUE)

auc_corrected   <- auc_apparent       - mean_opt_auc
slope_corrected <- cal_slope_apparent - mean_opt_slope
auc_m0_corrected <- auc_m0_apparent   - mean_opt_auc_m0

cat("=============================================================\n")
cat(" OPTIMISM-CORRECTED PERFORMANCE\n")
cat("=============================================================\n")
cat(sprintf(" %-30s  %8s  %12s  %12s\n",
            "Metric", "Apparent", "Mean optimism", "Corrected"))
cat(sprintf(" %-30s  %8.4f  %12.4f  %12.4f\n",
            "AUC", auc_apparent, mean_opt_auc, auc_corrected))
cat(sprintf(" %-30s  %8.4f  %12s  %12s\n",
            "Calibration intercept",
            cal_int_apparent, "— (see note)", "—"))
cat(sprintf(" %-30s  %8.4f  %12.4f  %12.4f\n",
            "Calibration slope",
            cal_slope_apparent, mean_opt_slope, slope_corrected))
cat(sprintf(" AUC clinical model (m0) apparent  : %.4f\n", auc_m0_apparent))
cat(sprintf(" AUC clinical model (m0) corrected : %.4f\n", auc_m0_corrected))  
cat("\n")

# =============================================================================
# 7. SAVE MAIN OUTPUT TABLE
# =============================================================================

output_table <- data.frame(
  Metric = c("AUC (clinical + SPAG1)",
             "AUC (clinical only — m0)",
             "Calibration intercept",
             "Calibration slope"),
  Ideal_value = c(NA_real_, NA_real_, 0, 1),
  Apparent = round(c(auc_apparent, auc_m0_apparent,
                     cal_int_apparent, cal_slope_apparent), 4),
  Mean_optimism = c(round(mean_opt_auc, 4),
                    round(mean_opt_auc_m0, 4),
                    NA_real_,
                    round(mean_opt_slope, 4)),
  Optimism_corrected = c(round(auc_corrected, 4),
                         round(auc_m0_corrected, 4),
                         NA_real_,
                         round(slope_corrected, 4)),
  Bootstrap_resamples = c(B, B, NA_integer_, B),
  Failed_resamples    = c(n_failed, n_failed, NA_integer_, n_failed),
  Seed                = c(SEED, SEED, NA_integer_, SEED),
  Note = c("", "", 
           "Algebraically 0 on training data; bootstrap not applicable",
           ""),
  stringsAsFactors = FALSE
)

# Hosmer-Lemeshow summary appended as attribute in a separate row block
hl_summary <- data.frame(
  Metric              = "Hosmer-Lemeshow test",
  Ideal_value         = NA_real_,
  Apparent            = NA_real_,
  Mean_optimism       = NA_real_,
  Optimism_corrected  = NA_real_,
  Bootstrap_resamples = NA_integer_,
  Failed_resamples    = NA_integer_,
  Seed                = NA_integer_,
  Note = sprintf("Chi2 = %.4f, df = %d, p = %.4f  (g = %d groups)",
                 hl_test$statistic, hl_test$parameter,
                 hl_test$p.value,   HL_G),
  stringsAsFactors = FALSE
)

output_table_full <- rbind(output_table, hl_summary)

write.csv(
  output_table_full,
  "tables/12_validation/Table_internal_validation.csv",
  row.names = FALSE,
  na = ""
)

cat(" Output table saved: tables/12_validation/Table_internal_validation.csv\n\n")
cat(" Contents:\n")
print(output_table_full, row.names = FALSE)
cat("\n")

# =============================================================================
# 8. CALIBRATION PLOT  (apparent, decile method)
# =============================================================================

dat_model$pred_prob <- prob_final

# Shrink predictions by corrected slope for corrected calibration plot
prob_corrected <- plogis(lp_final * slope_corrected)
dat_model$pred_prob_corrected <- prob_corrected

breaks <- unique(quantile(dat_model$pred_prob,
                          probs = seq(0, 1, 0.1),
                          na.rm = TRUE))
breaks[1]              <- breaks[1]              - 1e-10
breaks[length(breaks)] <- breaks[length(breaks)] + 1e-10

dat_model$decile <- cut(dat_model$pred_prob,
                        breaks         = breaks,
                        include.lowest = TRUE,
                        right          = TRUE)

cal_plot_df <- dat_model %>%
  group_by(decile) %>%
  summarise(
    mean_pred = mean(pred_prob),
    obs_rate  = mean(N_binary),
    n         = n(),
    .groups   = "drop"
  )

# HL p-value annotation for plot
hl_label <- sprintf("Hosmer-Lemeshow p = %.3f", hl_test$p.value)

p_cal <- ggplot(cal_plot_df, aes(x = mean_pred, y = obs_rate)) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", colour = "grey50", linewidth = 0.6) +
  geom_line(colour = "#2166AC", linewidth = 0.8) +
  geom_point(aes(size = n), colour = "#2166AC", fill = "#2166AC",
             shape = 21, alpha = 0.85) +
  # HL p-value annotation in top-left corner
  annotate("text",
           x = 0.02, y = 0.95,
           label = hl_label,
           hjust = 0, vjust = 1,
           size = 3.5, colour = "grey30") +
  scale_size_continuous(
    name   = "n per decile",
    range  = c(2, 7),
    breaks = c(30, 42, 60)
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = 0.02)
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = 0.02)
  ) +
  coord_equal() +
  labs(
    x       = "Mean predicted probability of N1",
    y       = "Observed N1 proportion",
    title   = "Calibration plot — apparent (decile method)",
    caption = sprintf(
      "Model: T-stage + Gleason group + SPAG1  |  n = %d  |  N1 events = %d (%.1f%%)",
      nrow(dat_model),
      sum(dat_model$N_binary),
      mean(dat_model$N_binary) * 100
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.caption     = element_text(colour = "grey40", size = 9),
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )

ggsave(
  "figures/12/Figure_calibration_plot_main_model.png",
  p_cal, width = 6, height = 6, dpi = 300
)

cat(" Calibration plot saved: figures/12/Figure_calibration_plot_main_model.png\n")

# ---- Corrected calibration plot ----
cal_plot_df_corr <- dat_model %>%
  mutate(decile_corr = cut(
    pred_prob_corrected,
    breaks = unique(c(
      min(pred_prob_corrected) - 1e-10,
      quantile(pred_prob_corrected, probs = seq(0.1, 0.9, 0.1)),
      max(pred_prob_corrected) + 1e-10
    )),
    include.lowest = TRUE
  )) %>%
  group_by(decile_corr) %>%
  summarise(
    mean_pred = mean(pred_prob_corrected),
    obs_rate  = mean(N_binary),
    n         = n(),
    .groups   = "drop"
  )

p_cal_corr <- ggplot(cal_plot_df_corr, aes(x = mean_pred, y = obs_rate)) +
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", colour = "grey50", linewidth = 0.6) +
  geom_line(colour = "#2166AC", linewidth = 0.8) +
  geom_point(aes(size = n), colour = "#2166AC", fill = "#2166AC",
             shape = 21, alpha = 0.85) +
  scale_size_continuous(
    name   = "n per decile",
    range  = c(2, 7),
    breaks = c(30, 42, 60)
  ) +
  scale_x_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = 0.02)
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = 0.02)
  ) +
  coord_equal() +
  labs(
    x       = "Mean corrected predicted probability of N1",
    y       = "Observed N1 proportion",
    title   = "Calibration plot — optimism-corrected (decile method)",
    caption = sprintf(
      "Model: T-stage + Gleason group + SPAG1  |  Calibration slope correction = %.4f  |  n = %d",
      slope_corrected, nrow(dat_model)
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title       = element_text(face = "bold", size = 13),
    plot.caption     = element_text(colour = "grey40", size = 9),
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )

ggsave(
  "figures/12/Figure_calibration_plot_corrected.png",
  p_cal_corr, width = 6, height = 6, dpi = 300
)
cat(" Corrected calibration plot saved: figures/12/Figure_calibration_plot_corrected.png\n")

# =============================================================================
# 9. COMBINED CALIBRATION FIGURE
# =============================================================================


p_cal_A <- p_cal +
  labs(
    title = "A. Apparent calibration",
    x = "Mean predicted probability of N1",
    y = "Observed N1 proportion"
  ) +
  theme(legend.position = "none")

p_cal_corr_B <- p_cal_corr +
  labs(
    title = "B. Optimism-corrected calibration",
    x = "Mean corrected predicted probability of N1",
    y = "Observed N1 proportion"
  )

combined_calib <- p_cal_A | p_cal_corr_B

ggsave(
  "figures/12/Figure_5_calibration_combined.png",
  combined_calib,
  width = 12,
  height = 6,
  dpi = 300
)

cat(" Combined calibration figure saved: figures/12/Figure_5_calibration_combined.png\n")
cat("=============================================================\n")
cat(" Script 12 complete.\n")
cat("=============================================================\n")

