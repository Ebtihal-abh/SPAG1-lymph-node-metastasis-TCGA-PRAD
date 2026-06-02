# ============================================================
# 08c_TSS_balance_SPAG1_quartiles.R
#
# Purpose:
#   Test whether tissue source site (TSS) is distributed
#   differently across SPAG1 Q1 vs Q4 groups.
#
#   Documents that TSS does not confound SPAG1 quartile
#   grouping, supporting the validity of the unadjusted
#   primary GSEA analysis in Script 08.
#
#   Mirrors the TSS balance check performed in Script 08b
#   for the N1 vs N0 analysis.
#
# Input:
#   data_processed/prad_master_features.rds 
#
# Output:
#   tables/08c/TSS_SPAG1_quartile_balance.csv
#   tables/08c/TSS_balance_summary.txt
# ============================================================

source("r/00_utils.R")  # loads collapse_tss() and TSS_MIN_N

library(dplyr)

# ── 1. Paths ─────────────────────────────────────────────────
data_path <- "data_processed/prad_master_features.rds"
out_dir   <- "tables/08c/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── 2. Load data ─────────────────────────────────────────────
dat <- readRDS(data_path)
cat("Master dataset loaded:", nrow(dat), "patients\n")
stopifnot("SPAG1"               %in% names(dat))
stopifnot("tissue_source_site"  %in% names(dat))

# ── 3. Define Q1 / Q4 groups (identical to Script 08) ────────
q1 <- quantile(dat$SPAG1, 0.25, na.rm = TRUE)
q3 <- quantile(dat$SPAG1, 0.75, na.rm = TRUE)

dat <- dat %>%
  mutate(
    SPAG1_group_q = case_when(
      SPAG1 <= q1 ~ "Low_Q1",
      SPAG1 >= q3 ~ "High_Q4",
      TRUE        ~ NA_character_
    )
  )

dat_sub <- dat %>% filter(!is.na(SPAG1_group_q))
cat("Q1 + Q4 patients retained:", nrow(dat_sub), "\n")
print(table(dat_sub$SPAG1_group_q))

# ── 4. Collapse TSS (same threshold as Scripts 08b and 13) ───
dat_sub <- dat_sub %>%
  mutate(TSS_collapsed = collapse_tss(tissue_source_site))

cat("\nTSS levels after collapsing (n >=", TSS_MIN_N, "):",
    nlevels(dat_sub$TSS_collapsed), "\n")
print(table(dat_sub$TSS_collapsed))

# ── 5. Cross-tabulate TSS vs SPAG1 quartile group ────────────
tss_tab <- table(
  TSS      = dat_sub$TSS_collapsed,
  SPAG1_Q  = dat_sub$SPAG1_group_q
)

cat("\nContingency table — TSS vs SPAG1 quartile group:\n")
print(tss_tab)

# ── 6. Chi-square test ────────────────────────────────────────
set.seed(123)
chi_res <- chisq.test(tss_tab, simulate.p.value = TRUE, B = 2000)

cat("\nChi-square test (TSS vs SPAG1 Q1/Q4):\n")
cat("  X-squared =", round(chi_res$statistic, 3), "\n")
cat("  p-value   =", signif(chi_res$p.value, 3),
    "(simulated, B = 2000)\n")

# ── 7. Per-TSS N1 rate across SPAG1 groups ───────────────────
tss_balance <- dat_sub %>%
  group_by(TSS_collapsed, SPAG1_group_q) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = SPAG1_group_q,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    total      = Low_Q1 + High_Q4,
    pct_High   = round(100 * High_Q4 / total, 1)
  ) %>%
  arrange(desc(total))

cat("\nPer-TSS distribution across SPAG1 Q1/Q4 groups:\n")
print(tss_balance)

write.csv(tss_balance,
          paste0(out_dir, "TSS_SPAG1_quartile_balance.csv"),
          row.names = FALSE)
cat("\nSaved:", paste0(out_dir, "TSS_SPAG1_quartile_balance.csv"), "\n")

# ── 8. Summary ───────────────────────────────────────────────
summary_lines <- c(
  "TSS Balance Check — SPAG1 Q1 vs Q4",
  paste0("Date: ", Sys.Date()),
  "",
  paste0("Q1 cutoff (25th percentile): ", round(q1, 4)),
  paste0("Q4 cutoff (75th percentile): ", round(q3, 4)),
  paste0("Patients in Q1 + Q4:         ", nrow(dat_sub)),
  paste0("TSS levels after collapsing: ", nlevels(dat_sub$TSS_collapsed)),
  "",
  paste0("Chi-square X-squared: ", round(chi_res$statistic, 3)),
  paste0("Chi-square p-value:   ", signif(chi_res$p.value, 3),
         " (simulated, B = 2000)"),
  "",
  paste0("Conclusion: TSS is ",
         ifelse(chi_res$p.value < 0.05, "", "NOT "),
         "associated with SPAG1 quartile grouping. ",
         ifelse(chi_res$p.value >= 0.05,
                "Unadjusted GSEA is not confounded by tissue source site.",
                "Consider TSS adjustment in GSEA sensitivity analysis."))
)

writeLines(summary_lines,
           paste0(out_dir, "TSS_balance_summary.txt"))
cat("Saved:", paste0(out_dir, "TSS_balance_summary.txt"), "\n")
cat("\n════ Script 08c Complete ════\n")