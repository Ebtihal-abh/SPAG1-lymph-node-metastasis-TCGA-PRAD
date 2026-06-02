# ============================================================
# Script 15: AUC Comparison — Clinical vs Clinical + SPAG1
# ============================================================
# Study:   TCGA-PRAD RNA-seq cohort (Illumina HiSeqV2)
# Outcome: Lymph node involvement (N1 = 1, N0 = 0)
# Purpose: Compare discrimination of clinical-only model (m0)
#          vs clinical + SPAG1 model (m1) using paired DeLong
#          test and bootstrap percentile CI for ΔAUC.
#
# Inputs:
#   data_processed/prad_master_features.rds
#
# Outputs:
#   tables/15/AUC_models.csv
#   tables/15/AUC_delta_and_pvalue.csv
#   tables/15/AUC_delta_bootstrap_CI.csv
#   figures/15/Figure_ROC_comparison.png
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(pROC)
  library(ggplot2)
  library(scales)
})

dir.create("tables/15", recursive = TRUE, showWarnings = FALSE)

# ---- Load master ----
dat <- readRDS("data_processed/prad_master_features.rds")

# ---- N_binary ----
dat <- dat %>%
  mutate(
    N_binary = case_when(
      tolower(as.character(N_group2)) == "n1" ~ 1L,
      tolower(as.character(N_group2)) == "n0" ~ 0L,
      TRUE ~ NA_integer_
    )
  )

# ---- Define variables (must match modeling set) ----
vars <- c("N_binary", "T_group2", "Gleason_group2", "SPAG1")

# ---- Complete-case dataset ----
dat_cc <- dat[complete.cases(dat[, vars]), , drop = FALSE]

cat("Samples used:", nrow(dat_cc), "\n")
cat("Events N1:", sum(dat_cc$N_binary == 1), " | N0:", sum(dat_cc$N_binary == 0), "\n")

# Ensure predictors are treated correctly
dat_cc <- dat_cc %>%
  mutate(
    T_group2 = factor(T_group2),
    Gleason_group2 = factor(Gleason_group2)
  )

# ---- Fit models ----
m0 <- glm(N_binary ~ T_group2 + Gleason_group2,
          data = dat_cc, family = binomial)

m1 <- glm(N_binary ~ T_group2 + Gleason_group2 + SPAG1,
          data = dat_cc, family = binomial)

# ---- Predicted probabilities ----
p0 <- predict(m0, type = "response")
p1 <- predict(m1, type = "response")

# ---- ROC + AUC (paired) ----
roc0 <- roc(response = dat_cc$N_binary, predictor = p0, quiet = TRUE)
roc1 <- roc(response = dat_cc$N_binary, predictor = p1, quiet = TRUE)

auc0 <- as.numeric(auc(roc0))
auc1 <- as.numeric(auc(roc1))
delta_auc <- auc1 - auc0

# DeLong test for paired ROCs (same patients, different models)
delong <- roc.test(roc0, roc1, method = "delong", paired = TRUE)

# ---- Output table (main result) ----
out_tbl <- data.frame(
  Model = c("Clinical_only (T + Gleason)", "Clinical + SPAG1 (T + Gleason + SPAG1)"),
  AUC   = c(auc0, auc1),
  stringsAsFactors = FALSE
)

delta_tbl <- data.frame(
  Delta_AUC = delta_auc,
  DeLong_p_value = as.numeric(delong$p.value),
  stringsAsFactors = FALSE
)

write.csv(out_tbl,   "tables/15/AUC_models.csv", row.names = FALSE)
write.csv(delta_tbl, "tables/15/AUC_delta_and_pvalue.csv", row.names = FALSE)

cat("\nAUC clinical-only:", round(auc0, 4), "\n")
cat("AUC + SPAG1      :", round(auc1, 4), "\n")
cat("Delta AUC        :", round(delta_auc, 4), "\n")
cat("DeLong p-value   :", format.pval(delong$p.value, digits = 3), "\n")

# bootstrap CI for Delta AUC ----
set.seed(123)
B <- 2000  
n <- nrow(dat_cc)

boot_delta <- replicate(B, {
  idx <- sample.int(n, size = n, replace = TRUE)
  d <- dat_cc[idx, , drop = FALSE]
  
  m0b <- glm(N_binary ~ T_group2 + Gleason_group2, data = d, family = binomial)
  m1b <- glm(N_binary ~ T_group2 + Gleason_group2 + SPAG1, data = d, family = binomial)
  
  p0b <- predict(m0b, type = "response")
  p1b <- predict(m1b, type = "response")
  
  r0 <- roc(d$N_binary, p0b, quiet = TRUE)
  r1 <- roc(d$N_binary, p1b, quiet = TRUE)
  
  as.numeric(auc(r1) - auc(r0))
})

ci <- quantile(boot_delta, probs = c(0.025, 0.975), na.rm = TRUE)

boot_tbl <- data.frame(
  Delta_AUC = delta_auc,
  CI_low = as.numeric(ci[1]),
  CI_high = as.numeric(ci[2]),
  B = B,
  stringsAsFactors = FALSE
)

write.csv(boot_tbl, "tables/15/AUC_delta_bootstrap_CI.csv", row.names = FALSE)

cat("Bootstrap 95% CI for ΔAUC:", round(ci[1], 4), "to", round(ci[2], 4), "\n")
cat("\nSaved tables in tables/15/\n")

# ============================================================
# ROC comparison figure
# ============================================================
dir.create("figures/15", recursive = TRUE, showWarnings = FALSE)

# Build ROC curve data frames from existing roc0 and roc2 objects
roc_df0 <- data.frame(
  fpr       = 1 - roc0$specificities,
  tpr       = roc0$sensitivities,
  Model     = sprintf("Clinical only  (AUC = %.3f)", auc0)
)

roc_df1 <- data.frame(
  fpr       = 1 - roc1$specificities,
  tpr       = roc1$sensitivities,
  Model     = sprintf("Clinical + SPAG1  (AUC = %.3f)", auc1)
)

roc_plot_df <- bind_rows(roc_df0, roc_df1) %>%
  mutate(Model = factor(Model, levels = c(
    sprintf("Clinical only  (AUC = %.3f)", auc0),
    sprintf("Clinical + SPAG1  (AUC = %.3f)", auc1)
  )))

# Annotation text
annot_label <- sprintf(
  "\u0394AUC = %.3f (95%% CI %.3f\u2013%.3f)\nDeLong p = 3.03\u00d710\u207b\u2075",
  delta_auc,
  as.numeric(ci[1]),
  as.numeric(ci[2])
)

# Build colour mapping before ggplot call
model_colours <- c("#999999", "#2166AC")
names(model_colours) <- c(
  sprintf("Clinical only  (AUC = %.3f)", auc0),
  sprintf("Clinical + SPAG1  (AUC = %.3f)", auc1)
)

p_roc <- ggplot(roc_plot_df,
                aes(x = fpr, y = tpr, colour = Model)) +
  geom_abline(intercept = 0, slope = 1,
              linetype  = "dashed",
              colour    = "grey60",
              linewidth = 0.5) +
  geom_line(linewidth = 0.9) +
  annotate("text",
           x      = 0.45,
           y      = 0.25,
           label  = annot_label,
           size   = 3.4,
           colour = "grey25",
           hjust  = 0) +
  scale_colour_manual(values = model_colours) +
  scale_x_continuous(
    name   = "1 \u2212 Specificity",
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = 0.02)
  ) +
  scale_y_continuous(
    name   = "Sensitivity",
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = expansion(mult = 0.02)
  ) +
  coord_equal() +
  labs(
    title   = "ROC curves \u2014 nodal metastasis prediction",
    colour  = NULL,
    caption = sprintf(
      "Model: T-stage + Gleason grade \u00b1 SPAG1  |  n = %d  |  N1 events = %d",
      nrow(dat_cc), sum(dat_cc$N_binary == 1)
    )
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position   = c(0.68, 0.15),
    legend.background = element_rect(fill      = "white",
                                     colour    = "grey80",
                                     linewidth = 0.4),
    legend.text       = element_text(size = 10),
    legend.key.width  = unit(1.2, "cm"),
    plot.title        = element_text(face = "bold", size = 13),
    plot.caption      = element_text(colour = "grey40", size = 9),
    panel.grid.minor  = element_blank()
  )

ggsave(
  "figures/15/Figure_ROC_comparison.png",
  p_roc,
  width  = 6,
  height = 6,
  dpi    = 300
)
message("Saved: figures/15/Figure_ROC_comparison.png")

cat("=============================================================\n")
cat(" Script 15 complete.\n")
cat("=============================================================\n")