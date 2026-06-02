# ============================================
# 07_biology_figures.R
# ============================================
# Generates Supplementary Figure S1 & Figure S2 plus a tidy
# CSV of panel statistics
# ============================================

library(dplyr)
library(ggplot2)
library(patchwork)

dir.create("figures/07", showWarnings = FALSE, recursive = TRUE)
dir.create("tables/07",  showWarnings = FALSE, recursive = TRUE)

dat <- readRDS("data_processed/prad_master_features.rds")

# ---- Helper: compute summary + test statistic ----

stat_two_group <- function(values, group, panel_label, var_label) {
  d <- data.frame(values = values, group = as.character(group)) %>%
    filter(!is.na(values), !is.na(group), group != "")
  
  groups <- sort(unique(d$group))
  if (length(groups) != 2) stop("stat_two_group expects exactly 2 levels")
  
  v1 <- d$values[d$group == groups[1]]
  v2 <- d$values[d$group == groups[2]]
  
  test <- suppressWarnings(wilcox.test(v1, v2))
  
  data.frame(
    panel = panel_label,
    variable = var_label,
    test = "Wilcoxon rank-sum",
    n_group1 = length(v1),
    median_group1 = round(median(v1), 2),
    Q1_group1 = round(quantile(v1, 0.25), 2),
    Q3_group1 = round(quantile(v1, 0.75), 2),
    label_group1 = groups[1],
    n_group2 = length(v2),
    median_group2 = round(median(v2), 2),
    Q1_group2 = round(quantile(v2, 0.25), 2),
    Q3_group2 = round(quantile(v2, 0.75), 2),
    label_group2 = groups[2],
    statistic = unname(test$statistic),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
}

stat_correlation <- function(x, y, panel_label, var_label) {
  d <- data.frame(x = x, y = y) %>% filter(!is.na(x), !is.na(y))
  ct <- suppressWarnings(cor.test(d$x, d$y, method = "spearman"))
  
  data.frame(
    panel = panel_label,
    variable = var_label,
    test = "Spearman correlation",
    n = nrow(d),
    rho = round(unname(ct$estimate), 3),
    statistic = NA_real_,
    p_value = ct$p.value,
    stringsAsFactors = FALSE
  )
}

# ---- Helper: format p-value as plotmath expression for figure annotation ----
format_p_label <- function(p) {
  if (p < 0.001) {
    exp_part  <- floor(log10(p))
    coef_part <- p / 10^exp_part
    sprintf("italic(p) == %.1f %%*%% 10^%d", coef_part, exp_part)
  } else {
    sprintf("italic(p) == %.3f", p)
  }
}

format_corr_label <- function(rho, p) {
  if (p < 0.001) {
    exp_part  <- floor(log10(p))
    coef_part <- p / 10^exp_part
    sprintf("rho == %.2f * ',' ~ italic(p) == %.1f %%*%% 10^%d",
            rho, coef_part, exp_part)
  } else {
    sprintf("rho == %.2f * ',' ~ italic(p) == %.3f", rho, p)
  }
}

# ---- Compute statistics for each panel ----

stat_A <- stat_two_group(dat$SPAG1, dat$N_group2,
                         "A", "SPAG1 by N stage")
stat_B <- stat_two_group(dat$SPAG1, dat$T_group2,
                         "B", "SPAG1 by T stage")
stat_C <- stat_two_group(dat$SPAG1, dat$Gleason_group2,
                         "C", "SPAG1 by Gleason group")
stat_D <- stat_correlation(dat$ImmuneScore, dat$SPAG1,
                           "A", "SPAG1 vs ImmuneScore")
stat_E <- stat_correlation(dat$StromalScore, dat$SPAG1,
                           "B", "SPAG1 vs StromalScore")

# TIMER2.0 purity result — pre-computed from web portal
# https://timer.cistrome.org/ — accessed April 2026
# Raw purity estimates not available for plotting
timer_rho <- -0.017
timer_p   <- 0.711

# Combine the two-group stats into one tidy block
two_group <- bind_rows(stat_A, stat_B, stat_C)
timer_row <- data.frame(
  panel     = "F",
  variable  = "SPAG1 vs Tumour Purity (TIMER2.0)",
  test      = "Spearman correlation (pre-computed)",
  n         = 497,
  rho       = timer_rho,
  statistic = NA_real_,
  p_value   = timer_p,
  stringsAsFactors = FALSE
)
correl <- bind_rows(stat_D, stat_E, timer_row)

# Save separately because column structures differ
write.csv(two_group, "tables/07/Fig_S1_panel_statistics_groupwise.csv",
          row.names = FALSE)
write.csv(correl,    "tables/07/Fig_S1_panel_statistics_correlation.csv",
          row.names = FALSE)

cat("\n=== Panel statistics for Supp Fig S1 ===\n\n")
cat("--- Panels A, B, C (Wilcoxon rank-sum) ---\n")
print(two_group, row.names = FALSE)
cat("\n--- Panels A, B (Spearman correlation) ---\n")
print(correl, row.names = FALSE)
cat("\n")

# ============================================
# 1) SPAG1 by N stage (Panel A)
# ============================================
p4 <- dat %>% filter(!is.na(N_group2)) %>%
  ggplot(aes(x = N_group2, y = SPAG1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  annotate("text", x = 1.5, y = Inf, vjust = 1.5,                                     
           label = format_p_label(stat_A$p_value),                                    
           parse = TRUE, size = 4) +                                                  
  labs(x = "Pathologic N group",
       y = "SPAG1 expression",
       title = "SPAG1 by N stage")

# ============================================
# 2) SPAG1 by T stage (Panel B)
# ============================================
p3 <- dat %>% filter(!is.na(T_group2)) %>%
  ggplot(aes(x = T_group2, y = SPAG1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  annotate("text", x = 1.5, y = Inf, vjust = 1.5,                                     
           label = format_p_label(stat_B$p_value),                                    
           parse = TRUE, size = 4) +                                                  
  labs(x = "Pathologic T group",
       y = "SPAG1 expression",
       title = "SPAG1 by T stage")



# ============================================
# 3) SPAG1 by Gleason group (Panel C)
# ============================================
p3b <- dat %>% filter(!is.na(Gleason_group2)) %>%
  ggplot(aes(x = Gleason_group2, y = SPAG1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, alpha = 0.4) +
  annotate("text", x = 1.5, y = Inf, vjust = 1.5,                                     
           label = format_p_label(stat_C$p_value),                                    
           parse = TRUE, size = 4) +                                                  
  labs(x = "Gleason group",
       y = "SPAG1 expression",
       title = "SPAG1 by Gleason grade")
# ============================================
# 4) SPAG1 vs ImmuneScore (Panel A)
# ============================================
p4b <- ggplot(dat, aes(x = ImmuneScore, y = SPAG1)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "steelblue") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = format_corr_label(stat_D$rho, stat_D$p_value),
           parse = TRUE, size = 4) +
  labs(x = "ESTIMATE ImmuneScore",
       y = "SPAG1 expression",
       title = "SPAG1 vs ImmuneScore")

# ============================================
# 5) SPAG1 vs StromalScore (Panel B)
# ============================================
p5 <- ggplot(dat, aes(x = StromalScore, y = SPAG1)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE, color = "darkgreen") +
  annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
           label = format_corr_label(stat_E$rho, stat_E$p_value),
           parse = TRUE, size = 4) +
  labs(x = "ESTIMATE StromalScore",
       y = "SPAG1 expression",
       title = "SPAG1 vs StromalScore")


# ============================================
# Figure 1: Clinical patterns (Panels A, B, C)
# ============================================
p4  <- p4  + ggtitle("A. SPAG1 by N stage")
p3  <- p3  + ggtitle("B. SPAG1 by T stage")
p3b <- p3b + ggtitle("C. SPAG1 by Gleason grade")

fig_clinical <- (p4 | p3 | p3b)

ggsave("figures/07/Fig_S1_clinical_patterns.png",
       fig_clinical,
       width = 14, height = 5, dpi = 300)
cat("Saved: figures/07/Fig_S1_clinical_patterns.png\n")

# ============================================
# Figure 2: Microenvironment correlations (Panels A, B)
# ============================================
p4b <- p4b + ggtitle("A. SPAG1 vs ImmuneScore")
p5  <- p5  + ggtitle("B. SPAG1 vs StromalScore")

fig_micro <- (p4b | p5)

ggsave("figures/07/Fig_S2_microenvironment.png",
       fig_micro,
       width = 10, height = 5, dpi = 300)

cat("\nAll figures and statistics saved in figures/07/ and tables/07/\n")