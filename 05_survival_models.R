# ============================================
# 05_survival_models.R
# Purpose:
# - Cox models for PFI
# - Save publication-ready HR tables
# - Check proportional hazards
# - KM curve for SPAG1 (median split for visualization only)
# - Martingale residual plot for linearity check (diagnostic only)
# ============================================

library(dplyr)
library(survival)
library(broom)
library(survminer)

dir.create("tables/05", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/05", showWarnings = FALSE)

# ---- Load dataset with features ----
dat <- readRDS("data_processed/prad_master_features.rds")

# ============================================================
# Survival endpoint event summary — justification for PFI choice
# ============================================================

endpoint_summary <- data.frame(
  Endpoint = c("OS", "DSS", "PFI", "DFI"),
  Time_var = c("OS.time", "DSS.time", "PFI.time", "DFI.time"),
  Event_var = c("OS", "DSS", "PFI", "DFI")
)

event_tbl <- lapply(1:nrow(endpoint_summary), function(i) {
  t_var <- endpoint_summary$Time_var[i]
  e_var <- endpoint_summary$Event_var[i]
  
  # Only patients with non-missing time AND event
  ok <- !is.na(dat[[t_var]]) & !is.na(dat[[e_var]])
  n_total  <- sum(ok)
  n_events <- sum(dat[[e_var]][ok] == 1, na.rm = TRUE)
  pct      <- round(100 * n_events / n_total, 1)
  
  data.frame(
    Endpoint     = endpoint_summary$Endpoint[i],
    N_evaluable  = n_total,
    N_events     = n_events,
    Event_pct    = pct,
    stringsAsFactors = FALSE
  )
})

event_tbl <- do.call(rbind, event_tbl)
event_tbl <- event_tbl[order(-event_tbl$N_events), ]

print(event_tbl)
write.csv(event_tbl, "tables/05/Survival_endpoint_event_summary.csv", row.names = FALSE)

# ---- Create survival variables (Xena uses PFI and PFI.time) ----
dat <- dat %>%
  mutate(
    PFI_time  = as.numeric(PFI.time),
    PFI_event = as.integer(PFI)
  )

# Sanity checks
stopifnot(all(dat$PFI_event %in% c(0,1)))
stopifnot(!all(is.na(dat$PFI_time)))

# ── Follow-up summary ─────────────────────────────────────────
followup_summary <- dat %>%
  summarise(
    n               = n(),
    median_followup = round(median(PFI_time / 365.25, na.rm = TRUE), 2),
    IQR_low         = round(quantile(PFI_time / 365.25, 0.25, na.rm = TRUE), 2),
    IQR_high        = round(quantile(PFI_time / 365.25, 0.75, na.rm = TRUE), 2),
    max_followup    = round(max(PFI_time / 365.25, na.rm = TRUE), 2),
    n_beyond_10yr   = sum(PFI_time / 365.25 > 10, na.rm = TRUE)
  )

cat(sprintf("Median follow-up : %.2f years (IQR %.2f–%.2f)\n",
            followup_summary$median_followup,
            followup_summary$IQR_low,
            followup_summary$IQR_high))
cat(sprintf("Maximum follow-up: %.2f years\n", followup_summary$max_followup))
cat(sprintf("Patients beyond 10 years: %d\n\n", followup_summary$n_beyond_10yr))

write.csv(followup_summary,
          "tables/05/PFI_followup_summary.csv",
          row.names = FALSE)
cat("Saved: tables/05/PFI_followup_summary.csv\n\n")
# --------------------------------------------
# 1) Survival analysis dataset (PFI available)
# --------------------------------------------
dat_surv_base <- dat %>%
  filter(!is.na(PFI_time),
         !is.na(PFI_event),
         !is.na(SPAG1))

cat("N with PFI + SPAG1:", nrow(dat_surv_base), "\n")
cat("PFI events:", sum(dat_surv_base$PFI_event == 1), "\n\n")

# --------------------------------------------
# 2) Models (PFI)
# --------------------------------------------

# Unadjusted
cox_unadj <- coxph(Surv(PFI_time, PFI_event) ~ SPAG1,
                   data = dat_surv_base)

# Fully adjusted: Gleason + T + N
dat_surv_GTN <- dat_surv_base %>%
  filter(!is.na(Gleason_group2),
         !is.na(T_group2),
         !is.na(N_group2)) %>%
  droplevels()

# Ensure reference levels 
dat_surv_GTN$Gleason_group2 <- factor(dat_surv_GTN$Gleason_group2)
dat_surv_GTN$T_group2 <- factor(dat_surv_GTN$T_group2)
dat_surv_GTN$N_group2 <- factor(dat_surv_GTN$N_group2, levels = c("N0","N1"))

cox_adj_GTN <- coxph(Surv(PFI_time, PFI_event) ~ SPAG1 + Gleason_group2 + T_group2 + N_group2,
                     data = dat_surv_GTN)

cat("N for unadjusted model:", nrow(dat_surv_base), "\n")
cat("PFI events (unadjusted):", sum(dat_surv_base$PFI_event == 1), "\n")
cat("N for fully adjusted model (GTN):", nrow(dat_surv_GTN), "\n")
cat("PFI events (GTN):", sum(dat_surv_GTN$PFI_event == 1), "\n")

# --------------------------------------------
# 3) Save tidy results tables (HR + CI)
# --------------------------------------------

tidy_cox <- function(fit) {
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    select(term, estimate, conf.low, conf.high, p.value)
}

write.csv(tidy_cox(cox_unadj),   "tables/05/PFI_cox_unadjusted.csv",            row.names = FALSE)
write.csv(tidy_cox(cox_adj_GTN), "tables/05/PFI_cox_adjusted_Gleason_T_N.csv",  row.names = FALSE)

# --------------------------------------------
# 4) Proportional hazards checks — ALL models
# --------------------------------------------

zph_unadj <- cox.zph(cox_unadj)
zph_GTN   <- cox.zph(cox_adj_GTN)

# Save to tables/05/ (consistent with all other outputs from this script)
capture.output(print(zph_unadj), file = "tables/05/PFI_PH_test_unadjusted.txt")
capture.output(print(zph_GTN),   file = "tables/05/PFI_PH_test_adjusted_Gleason_T_N.txt")


# Print global p-values to console for quick review
cat("\n--- Proportional hazards: GLOBAL p-values ---\n")
cat(sprintf("  Unadjusted:              p = %.3f\n", zph_unadj$table["GLOBAL", "p"]))
cat(sprintf("  Adjusted (Gleason+T+N):  p = %.3f\n", zph_GTN$table["GLOBAL", "p"]))

# ============================================================
# 5) Linearity check — Martingale residuals
# ============================================================
# Purpose: confirm linear relationship between SPAG1 expression
# and log-hazard in the Cox model.
# Diagnostic check only — not reported in manuscript.
# A flat lowess curve around zero confirms linearity.
# ============================================================

mart_resid <- residuals(cox_unadj, type = "martingale")
spag1_vals <- dat_surv_base$SPAG1

png("figures/05/Figure_linearity_check_cox.png",
    width = 6, height = 5, units = "in", res = 300)

plot(spag1_vals, mart_resid,
     xlab = "SPAG1 expression",
     ylab = "Martingale residuals",
     main = "Linearity check \u2014 Cox model",
     pch  = 1,
     col  = "grey40",
     cex  = 0.7)
abline(h = 0, col = "red", lwd = 1.2)
lines(lowess(spag1_vals, mart_resid), col = "steelblue", lwd = 2)

dev.off()
cat("Saved: figures/05/Figure_linearity_check_cox.png\n")

# Console summary
lowess_fit  <- lowess(spag1_vals, mart_resid)
resid_range <- range(lowess_fit$y)
cat(sprintf(
  "\nLinearity check: lowess range = [%.3f, %.3f]\n",
  resid_range[1], resid_range[2]
))
cat("Interpretation: flat lowess near zero confirms linear assumption.\n\n")

# ============================================================
# KM CURVE — SPAG1 High vs Low (median split, for visualization only)
# Note: Cox models above use SPAG1 as continuous — median split
#       is for visualization purposes only, as stated in figure legend
# ============================================================


# Median split on SPAG1
spag1_median <- median(dat_surv_base$SPAG1, na.rm = TRUE)

dat_km <- dat_surv_base %>%
  mutate(
    SPAG1_group = factor(
      ifelse(SPAG1 >= spag1_median, "High", "Low"),
      levels = c("Low", "High")
    )
  )

cat(sprintf("SPAG1 median: %.4f\n", spag1_median))
cat("KM group sizes:\n")
print(table(dat_km$SPAG1_group))
cat("Events per group:\n")
print(table(dat_km$SPAG1_group, dat_km$PFI_event))

# Fit KM
km_fit <- survfit(
  Surv(PFI_time, PFI_event) ~ SPAG1_group,
  data = dat_km
)

# Log-rank test
lr_test <- survdiff(
  Surv(PFI_time, PFI_event) ~ SPAG1_group,
  data = dat_km
)

lr_p <- 1 - pchisq(lr_test$chisq, df = length(lr_test$n) - 1)
cat(sprintf("\nLog-rank p-value: %.4f\n", lr_p))

# ---- Plot ----
p_km <- ggsurvplot(
  km_fit,
  data                = dat_km,
  pval                = TRUE,
  pval.method         = TRUE,
  conf.int            = TRUE,
  risk.table          = TRUE,
  risk.table.height   = 0.28,
  xscale              = "d_m",        # convert days to months
  break.time.by       = 365.25,       # break every 12 months
  xlim                = c(0, 3650),   # truncate at 10 years
  xlab                = "Time (months)",
  ylab                = "Progression-Free Interval probability",
  legend.title        = "SPAG1 expression",
  legend.labs         = c("Low (\u2264 median)", "High (> median)"),
  palette             = c("#2166AC", "#B2182B"),
  title               = "Kaplan-Meier: PFI by SPAG1 expression (median split)",
  ggtheme             = theme_bw(base_size = 12),
  font.main           = c(13, "bold"),
  risk.table.fontsize = 3.5,
  tables.theme        = theme_cleantable(),
  caption = sprintf(
    "n = %d  |  Events = %d (%.1f%%)  |  SPAG1 median = %.2f\nNote: median split for visualization only; Cox models use SPAG1 as continuous variable",
    nrow(dat_km), sum(dat_km$PFI_event),
    mean(dat_km$PFI_event) * 100, spag1_median
  )
)

# ---- Save using png/dev.off  ----

png("figures/05/Figure_KM_PFI_SPAG1.png",
    width = 8, height = 7, units = "in", res = 300)
print(p_km)
dev.off()

cat("Saved: figures/05/Figure_KM_PFI_SPAG1.png\n")

# Save log-rank result
write.csv(
  data.frame(
    Test        = "Log-rank",
    Endpoint    = "PFI",
    Group_var   = "SPAG1_median_split",
    SPAG1_median = round(spag1_median, 4),
    Chi_sq      = round(lr_test$chisq, 4),
    df          = length(lr_test$n) - 1,
    p_value     = round(lr_p, 4)
  ),
  "tables/05/KM_logrank_PFI_SPAG1.csv",
  row.names = FALSE
)
cat("Saved: tables/05/KM_logrank_PFI_SPAG1.csv\n")

