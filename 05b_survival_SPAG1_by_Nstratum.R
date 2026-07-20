# ============================================================
# 05b_survival_SPAG1_by_Nstratum.R
# Purpose:
# - Address reviewer request: association between SPAG1 and PFI
#   WITHIN node-positive (N1) tumours, with the node-negative (N0)
#   stratum reported alongside for symmetry.
# - Univariable Cox (SPAG1 continuous) within each nodal stratum.
#   NOTE: subgroups are small (esp. N1: ~79 patients, ~22 events),
#   so models are deliberately UNIVARIABLE, one predictor only,
#   consistent with events-per-variable constraints. No covariate
#   adjustment is attempted within strata.
# - KM curves (median split, visualization only) within each stratum.
# - Proportional-hazards checks for each model.
#
# Mirrors conventions of 05_survival_models.R.
#
# Input:
#   data_processed/prad_master_features.rds
#
# Outputs:
#   tables/05b/PFI_cox_SPAG1_within_N0.csv
#   tables/05b/PFI_cox_SPAG1_within_N1.csv
#   tables/05b/PFI_cox_SPAG1_by_Nstratum_summary.csv
#   tables/05b/PFI_PH_test_within_N0.txt
#   tables/05b/PFI_PH_test_within_N1.txt
#   tables/05b/KM_logrank_PFI_SPAG1_by_Nstratum.csv
#   figures/05b/Figure_KM_PFI_SPAG1_N0.png
#   figures/05b/Figure_KM_PFI_SPAG1_N1.png
# ============================================================

library(dplyr)
library(survival)
library(broom)
library(survminer)

dir.create("tables/05b",  recursive = TRUE, showWarnings = FALSE)
dir.create("figures/05b", recursive = TRUE, showWarnings = FALSE)

# ---- Load dataset with features ----
dat <- readRDS("data_processed/prad_master_features.rds")

# ---- Survival variables (Xena uses PFI and PFI.time) ----
dat <- dat %>%
  mutate(
    PFI_time  = as.numeric(PFI.time),
    PFI_event = as.integer(PFI),
    N_stratum = case_when(
      tolower(as.character(N_group2)) == "n1" ~ "N1",
      tolower(as.character(N_group2)) == "n0" ~ "N0",
      TRUE                                    ~ NA_character_
    )
  )

stopifnot(all(dat$PFI_event %in% c(0, 1) | is.na(dat$PFI_event)))

# ---- Base survival set: non-missing PFI, SPAG1, and nodal stratum ----
dat_surv <- dat %>%
  filter(!is.na(PFI_time), !is.na(PFI_event),
         !is.na(SPAG1), !is.na(N_stratum))

cat("=============================================================\n")
cat(" SPAG1 vs PFI — univariable Cox within nodal strata\n")
cat("=============================================================\n")

# ============================================================
# HELPER: univariable Cox (SPAG1 continuous) within one stratum
# ============================================================
tidy_cox <- function(fit) {
  broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
    select(term, estimate, conf.low, conf.high, p.value)
}

analyse_stratum_cox <- function(stratum) {
  
  d <- dat_surv %>% filter(N_stratum == stratum)
  n_pt <- nrow(d)
  n_ev <- sum(d$PFI_event == 1)
  
  cat("-------------------------------------------------------------\n")
  cat(sprintf(" Stratum %s: n = %d | PFI events = %d | EPV = %.1f\n",
              stratum, n_pt, n_ev, n_ev))  # EPV = events for 1 predictor
  
  # Univariable Cox: SPAG1 continuous (parallel to main models)
  fit <- coxph(Surv(PFI_time, PFI_event) ~ SPAG1, data = d)
  res <- tidy_cox(fit)
  
  write.csv(res,
            sprintf("tables/05b/PFI_cox_SPAG1_within_%s.csv", stratum),
            row.names = FALSE)
  
  # Proportional hazards check
  zph <- cox.zph(fit)
  capture.output(print(zph),
                 file = sprintf("tables/05b/PFI_PH_test_within_%s.txt", stratum))
  
  hr  <- res$estimate[res$term == "SPAG1"]
  lo  <- res$conf.low[res$term == "SPAG1"]
  hi  <- res$conf.high[res$term == "SPAG1"]
  p   <- res$p.value[res$term == "SPAG1"]
  cat(sprintf("   SPAG1 (continuous): HR = %.2f (95%% CI %.2f-%.2f), p = %.3f\n",
              hr, lo, hi, p))
  cat(sprintf("   PH test (Schoenfeld) p = %.3f\n",
              zph$table["SPAG1", "p"]))

  # ---- Linearity check — Martingale residuals (diagnostic) ----
  mart_resid <- residuals(fit, type = "martingale")
  png(sprintf("figures/05b/Figure_linearity_check_cox_%s.png", stratum),
      width = 6, height = 5, units = "in", res = 300)
  plot(d$SPAG1, mart_resid,
       xlab = "SPAG1 expression",
       ylab = "Martingale residuals",
       main = sprintf("Linearity check \u2014 Cox model (%s)", stratum),
       pch = 1, col = "grey40", cex = 0.7)
  abline(h = 0, col = "red", lwd = 1.2)
  lines(lowess(d$SPAG1, mart_resid), col = "steelblue", lwd = 2)
  dev.off()
  lin_range <- range(lowess(d$SPAG1, mart_resid)$y)
  cat(sprintf("   Linearity (lowess range): [%.3f, %.3f]\n",
              lin_range[1], lin_range[2]))
    
  # ---- KM (median split, visualization only) ----
  med <- median(d$SPAG1, na.rm = TRUE)
  d <- d %>%
    mutate(SPAG1_group = factor(ifelse(SPAG1 >= med, "High", "Low"),
                                levels = c("Low", "High")))
  
  km_fit  <- survfit(Surv(PFI_time, PFI_event) ~ SPAG1_group, data = d)
  lr_test <- survdiff(Surv(PFI_time, PFI_event) ~ SPAG1_group, data = d)
  lr_p    <- 1 - pchisq(lr_test$chisq, df = length(lr_test$n) - 1)
  
  p_km <- ggsurvplot(
    km_fit, data = d,
    pval = TRUE, conf.int = TRUE, risk.table = TRUE,
    risk.table.height = 0.28,
    xscale = "d_m", break.time.by = 365.25, xlim = c(0, 3650),
    xlab = "Time (months)", ylab = "Progression-Free Interval probability",
    legend.title = "SPAG1 expression",
    legend.labs  = c("Low (\u2264 median)", "High (> median)"),
    palette = c("#2166AC", "#B2182B"),
    title = sprintf("PFI by SPAG1 expression within %s (median split)", stratum),
    ggtheme = theme_bw(base_size = 12),
    risk.table.fontsize = 3.5, tables.theme = theme_cleantable(),
    caption = sprintf(
      "Stratum %s: n = %d | Events = %d | median split for visualization only",
      stratum, n_pt, n_ev)
  )
  
  png(sprintf("figures/05b/Figure_KM_PFI_SPAG1_%s.png", stratum),
      width = 8, height = 7, units = "in", res = 300)
  print(p_km)
  dev.off()
  
  data.frame(
    Stratum       = stratum,
    N             = n_pt,
    PFI_events    = n_ev,
    HR_SPAG1      = round(hr, 3),
    CI_low        = round(lo, 3),
    CI_high       = round(hi, 3),
    Cox_p         = signif(p, 3),
    PH_p          = signif(zph$table["SPAG1", "p"], 3),
    KM_logrank_p  = signif(lr_p, 3),
    SPAG1_median  = round(med, 4),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# Run both strata (N0 for symmetry, N1 = reviewer's question)
# ============================================================
summary_tbl <- bind_rows(lapply(c("N0", "N1"), analyse_stratum_cox))

write.csv(summary_tbl,
          "tables/05b/PFI_cox_SPAG1_by_Nstratum_summary.csv",
          row.names = FALSE)

# Log-rank summary (both strata)
write.csv(
  summary_tbl %>%
    select(Stratum, N, PFI_events, SPAG1_median, KM_logrank_p),
  "tables/05b/KM_logrank_PFI_SPAG1_by_Nstratum.csv",
  row.names = FALSE
)

cat("\n=============================================================\n")
cat(" SUMMARY — SPAG1 (continuous) vs PFI within strata\n")
cat("=============================================================\n")
print(summary_tbl, row.names = FALSE)

cat(sprintf(
  paste0("\n Note: the N1 subgroup is small (n = %d, %d events); the within-N1\n",
         " estimate is univariable and should be read as exploratory. A null or\n",
         " weak association is consistent with the whole-cohort adjusted model,\n",
         " in which SPAG1 was not an independent predictor of PFI.\n"),
  summary_tbl$N[summary_tbl$Stratum == "N1"],
  summary_tbl$PFI_events[summary_tbl$Stratum == "N1"]))

cat("\n\u2550\u2550\u2550\u2550 Script 05b complete \u2550\u2550\u2550\u2550\n")