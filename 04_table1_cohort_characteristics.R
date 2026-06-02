# ============================================
# 04_table1_cohort_characteristics.R
#
# Purpose:
#   Generate Table 1 (cohort characteristics by pathological N stage)
#
# Input:
#   prad_master_features.rds
#
# Output:
#   tables/Table1_cohort_characteristics.csv
# ============================================

library(dplyr)

# ---- Paths ----
path_in   <- "data_processed/prad_master_features.rds"
path_out  <- "tables/04/Table1_cohort_characteristics.csv"

dir.create(dirname(path_out), showWarnings = FALSE, recursive = TRUE)

# ---- 1) Load data ----
dat <- readRDS(path_in)
cat("Loaded:", nrow(dat), "patients,", ncol(dat), "variables\n\n")

# ---- 2) Define N stage grouping using N_group2 directly ----
# N_group2 has: N0 (345), N1 (79), NA (73 = NX)
dat$Nstage_table1 <- ifelse(
  is.na(dat$N_group2), "NX",
  as.character(dat$N_group2)
)

cat("N stage distribution for Table 1:\n")
print(table(dat$Nstage_table1, useNA = "ifany"))
cat("\n")

# ---- 3) Compute age in years ----
if ("days_to_birth" %in% names(dat)) {
  dat$age_years <- -as.numeric(dat$days_to_birth) / 365.25
} else if ("age_at_initial_pathologic_diagnosis" %in% names(dat)) {
  dat$age_years <- as.numeric(dat$age_at_initial_pathologic_diagnosis)
}

# ---- 4) Helper functions ----

fmt_p <- function(p) {
  if (is.na(p)) return("\u2014")
  if (p < 0.001) return(formatC(p, format = "e", digits = 1))
  return(formatC(p, format = "f", digits = 3))
}

fmt_iqr <- function(x, dp = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("\u2014")
  med <- formatC(median(x), format = "f", digits = dp)
  q1  <- formatC(quantile(x, 0.25, names = FALSE), format = "f", digits = dp)
  q3  <- formatC(quantile(x, 0.75, names = FALSE), format = "f", digits = dp)
  paste0(med, " (", q1, "\u2013", q3, ")")
}

fmt_n_pct <- function(n, total, dp = 1) {
  if (is.na(n) || total == 0) return("0 (0.0%)")
  pct <- 100 * n / total
  paste0(n, " (", formatC(pct, format = "f", digits = dp), "%)")
}

# Continuous test: Wilcoxon between N0 and N1
test_cont <- function(values, group) {
  keep <- group %in% c("N0", "N1") & !is.na(values)
  if (sum(keep) < 5) return(NA_real_)
  v0 <- values[keep & group == "N0"]
  v1 <- values[keep & group == "N1"]
  if (length(v0) < 2 || length(v1) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(v0, v1)$p.value)
}

# Categorical test: chi-square or Fisher (suppressed warnings when Fisher used)
test_cat <- function(values, group) {
  keep <- group %in% c("N0", "N1") & !is.na(values) & values != ""
  if (sum(keep) < 5) return(NA_real_)
  tab <- table(values[keep], group[keep])
  if (any(dim(tab) < 2)) return(NA_real_)
  expected <- suppressWarnings(chisq.test(tab)$expected)
  if (any(expected < 5)) {
    return(suppressWarnings(fisher.test(tab)$p.value))
  } else {
    return(suppressWarnings(chisq.test(tab)$p.value))
  }
}

# ---- 5) Build Table 1 ----

rows <- list()
add_row <- function(variable, n0, n1, nx, p) {
  rows[[length(rows) + 1]] <<- data.frame(
    Characteristic = variable,
    N0 = n0,
    N1 = n1,
    NX = nx,
    `p-value (N0 vs N1)` = p,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

g  <- dat$Nstage_table1
n0 <- sum(g == "N0", na.rm = TRUE)
n1 <- sum(g == "N1", na.rm = TRUE)
nx <- sum(g == "NX", na.rm = TRUE)

add_row("Total, n", as.character(n0), as.character(n1), as.character(nx), "")

# ---- Age ----
if ("age_years" %in% names(dat)) {
  p <- test_cont(dat$age_years, g)
  add_row(
    "Age at diagnosis, years, median (IQR)",
    fmt_iqr(dat$age_years[g == "N0"], dp = 0),
    fmt_iqr(dat$age_years[g == "N1"], dp = 0),
    fmt_iqr(dat$age_years[g == "NX"], dp = 0),
    fmt_p(p)
  )
}

# ---- PSA at sample collection ----
# These values appear to be post-operative (most are 0.05 ng/mL),
# not pre-operative diagnostic PSA. Labelled accordingly.
if ("psa_value" %in% names(dat)) {
  psa <- as.numeric(dat$psa_value)
  if (sum(!is.na(psa)) > nrow(dat) * 0.5) {
    p <- test_cont(psa, g)
    add_row(
      "PSA at sample collection, ng/mL, median (IQR)",
      fmt_iqr(psa[g == "N0"], dp = 2),
      fmt_iqr(psa[g == "N1"], dp = 2),
      fmt_iqr(psa[g == "NX"], dp = 2),
      fmt_p(p)
    )
  }
}

# ---- T stage ----
if ("T_group2" %in% names(dat)) {
  t_var <- as.character(dat$T_group2)
  p <- test_cat(t_var, g)
  add_row("Pathological T stage, n (%)", "", "", "", fmt_p(p))
  add_row("  T2",
          fmt_n_pct(sum(t_var == "T2"   & g == "N0", na.rm = TRUE), n0),
          fmt_n_pct(sum(t_var == "T2"   & g == "N1", na.rm = TRUE), n1),
          fmt_n_pct(sum(t_var == "T2"   & g == "NX", na.rm = TRUE), nx),
          "")
  add_row("  T3/T4",
          fmt_n_pct(sum(t_var == "T3_4" & g == "N0", na.rm = TRUE), n0),
          fmt_n_pct(sum(t_var == "T3_4" & g == "N1", na.rm = TRUE), n1),
          fmt_n_pct(sum(t_var == "T3_4" & g == "NX", na.rm = TRUE), nx),
          "")
  # Show missing T separately if any
  n_miss_t <- sum(is.na(t_var))
  if (n_miss_t > 0) {
    add_row("  Unknown",
            fmt_n_pct(sum(is.na(t_var) & g == "N0"), n0),
            fmt_n_pct(sum(is.na(t_var) & g == "N1"), n1),
            fmt_n_pct(sum(is.na(t_var) & g == "NX"), nx),
            "")
  }
}

# ---- Gleason group ----
if ("Gleason_group2" %in% names(dat)) {
  gl_var <- as.character(dat$Gleason_group2)
  p <- test_cat(gl_var, g)
  add_row("Gleason grade group, n (%)", "", "", "", fmt_p(p))
  add_row("  Gleason \u22647 (Low/Intermediate)",
          fmt_n_pct(sum(gl_var == "Low_intermediate" & g == "N0", na.rm = TRUE), n0),
          fmt_n_pct(sum(gl_var == "Low_intermediate" & g == "N1", na.rm = TRUE), n1),
          fmt_n_pct(sum(gl_var == "Low_intermediate" & g == "NX", na.rm = TRUE), nx),
          "")
  add_row("  Gleason \u22658 (High)",
          fmt_n_pct(sum(gl_var == "High"             & g == "N0", na.rm = TRUE), n0),
          fmt_n_pct(sum(gl_var == "High"             & g == "N1", na.rm = TRUE), n1),
          fmt_n_pct(sum(gl_var == "High"             & g == "NX", na.rm = TRUE), nx),
          "")
}

# ---- Surgical margin ----
if ("residual_tumor" %in% names(dat)) {
  rt <- as.character(dat$residual_tumor)
  rt[rt == ""] <- NA
  if (sum(!is.na(rt)) > nrow(dat) * 0.5) {
    rt_bin <- ifelse(rt == "R0", "Negative",
                     ifelse(rt %in% c("R1", "R2"), "Positive", NA))
    p <- test_cat(rt_bin, g)
    add_row("Surgical margin, n (%)", "", "", "", fmt_p(p))
    for (lvl in c("Negative", "Positive")) {
      add_row(paste0("  ", lvl),
              fmt_n_pct(sum(rt_bin == lvl & g == "N0", na.rm = TRUE), n0),
              fmt_n_pct(sum(rt_bin == lvl & g == "N1", na.rm = TRUE), n1),
              fmt_n_pct(sum(rt_bin == lvl & g == "NX", na.rm = TRUE), nx),
              "")
    }
  }
}

# ---- Biochemical recurrence ----
if ("biochemical_recurrence" %in% names(dat)) {
  bcr <- as.character(dat$biochemical_recurrence)
  bcr[bcr == ""] <- NA
  if (sum(!is.na(bcr)) > nrow(dat) * 0.5) {
    p <- test_cat(bcr, g)
    add_row("Biochemical recurrence, n (%)", "", "", "", fmt_p(p))
    add_row("  Yes",
            fmt_n_pct(sum(bcr == "YES" & g == "N0", na.rm = TRUE), n0),
            fmt_n_pct(sum(bcr == "YES" & g == "N1", na.rm = TRUE), n1),
            fmt_n_pct(sum(bcr == "YES" & g == "NX", na.rm = TRUE), nx),
            "")
    add_row("  No",
            fmt_n_pct(sum(bcr == "NO"  & g == "N0", na.rm = TRUE), n0),
            fmt_n_pct(sum(bcr == "NO"  & g == "N1", na.rm = TRUE), n1),
            fmt_n_pct(sum(bcr == "NO"  & g == "NX", na.rm = TRUE), nx),
            "")
  }
}

# ---- Vital status ----
if ("vital_status" %in% names(dat)) {
  vs <- as.character(dat$vital_status)
  vs[vs == ""] <- NA
  if (sum(!is.na(vs)) > nrow(dat) * 0.5) {
    p <- test_cat(vs, g)
    add_row("Vital status, n (%)", "", "", "", fmt_p(p))
    add_row("  Alive",
            fmt_n_pct(sum(vs == "LIVING"   & g == "N0", na.rm = TRUE), n0),
            fmt_n_pct(sum(vs == "LIVING"   & g == "N1", na.rm = TRUE), n1),
            fmt_n_pct(sum(vs == "LIVING"   & g == "NX", na.rm = TRUE), nx),
            "")
    add_row("  Deceased",
            fmt_n_pct(sum(vs == "DECEASED" & g == "N0", na.rm = TRUE), n0),
            fmt_n_pct(sum(vs == "DECEASED" & g == "N1", na.rm = TRUE), n1),
            fmt_n_pct(sum(vs == "DECEASED" & g == "NX", na.rm = TRUE), nx),
            "")
  }
}

# ---- Follow-up time (months) ----
fup_var <- NULL
fup_label <- NULL
if ("PFI.time" %in% names(dat)) {
  fup_var <- as.numeric(dat$PFI.time) / 30.4375
  fup_label <- "Follow-up time (PFI), months, median (IQR)"
} else if ("days_to_last_followup" %in% names(dat)) {
  fup_var <- as.numeric(dat$days_to_last_followup) / 30.4375
  fup_label <- "Follow-up time, months, median (IQR)"
}
if (!is.null(fup_var) && sum(!is.na(fup_var)) > nrow(dat) * 0.5) {
  p <- test_cont(fup_var, g)
  add_row(
    fup_label,
    fmt_iqr(fup_var[g == "N0"], dp = 1),
    fmt_iqr(fup_var[g == "N1"], dp = 1),
    fmt_iqr(fup_var[g == "NX"], dp = 1),
    fmt_p(p)
  )
}

# ---- PFI events ----
if ("PFI" %in% names(dat)) {
  pfi <- as.numeric(dat$PFI)
  add_row("Progression-free interval events, n (%)",
          fmt_n_pct(sum(pfi == 1 & g == "N0", na.rm = TRUE), n0),
          fmt_n_pct(sum(pfi == 1 & g == "N1", na.rm = TRUE), n1),
          fmt_n_pct(sum(pfi == 1 & g == "NX", na.rm = TRUE), nx),
          "")
}

# ---- SPAG1 expression ----
if ("SPAG1" %in% names(dat)) {
  p <- test_cont(dat$SPAG1, g)
  add_row(
    "SPAG1 expression (log\u2082 RSEM), median (IQR)",
    fmt_iqr(dat$SPAG1[g == "N0"], dp = 2),
    fmt_iqr(dat$SPAG1[g == "N1"], dp = 2),
    fmt_iqr(dat$SPAG1[g == "NX"], dp = 2),
    fmt_p(p)
  )
}

# ---- Combine and save ----
out <- do.call(rbind, rows)

cat("\n")
print(out, row.names = FALSE)

write.csv(out, path_out, row.names = FALSE, fileEncoding = "UTF-8")
cat("\nSaved Table 1 to:", path_out, "\n")
cat("Rows in Table 1:", nrow(out), "\n")