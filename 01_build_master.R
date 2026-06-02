library(dplyr)


# Load TCGA datasets
clin <- read.delim("data_raw/clinical-data.tsv", check.names = FALSE)
surv <- read.delim("data_raw/survival-data.tsv", check.names = FALSE)


# Clean, consistent TCGA patient barcode
# Retain Primary Tumor samples only before deduplication
clin2 <- clin %>%
  filter(sample_type == "Primary Tumor") %>%
  mutate(patient_id = substr(trimws(as.character(bcr_patient_barcode)), 1, 12)) %>%
  arrange(patient_id, bcr_patient_barcode) %>%        # deterministic row selection
  distinct(patient_id, .keep_all = TRUE)

cat("Clinical: rows before dedup:", sum(clin$sample_type == "Primary Tumor", na.rm=TRUE),
    "| after:", nrow(clin2),
    "| duplicates removed:", sum(clin$sample_type == "Primary Tumor", na.rm=TRUE) - nrow(clin2), "\n")

surv2 <- surv %>%
  mutate(patient_id = substr(trimws(as.character(sample)), 1, 12)) %>%
  arrange(patient_id, sample) %>%                     # deterministic row selection
  distinct(patient_id, .keep_all = TRUE)

cat("Survival: rows before dedup:", nrow(surv),
    "| after:", nrow(surv2),
    "| duplicates removed:", nrow(surv) - nrow(surv2), "\n\n")

# Sanity: how many match?
n_match <- length(intersect(clin2$patient_id, surv2$patient_id))
cat("Patients matched (clin + surv):", n_match, "\n")
cat("In clinical only:", nrow(clin2) - n_match, "\n")
cat("In survival only:", nrow(surv2) - n_match, "\n")

# Duplicate check results (verified):
# Clinical: 498 Primary Tumor samples, 0 duplicate patient_ids — no rows removed
# Survival: 566 rows (sample-level) → 498 unique patients after dedup
#           68 rows removed = patients with multiple sample entries in survival file
#          TCGA survival files are sample-level, not patient-level
# Match: 498/498 — complete overlap, no missing survival data after join

# Join
dat <- left_join(
  clin2,
  surv2 %>% select(patient_id, PFI, PFI.time, OS, OS.time, DSS, DSS.time, DFI, DFI.time),
  by = "patient_id"
)


dir.create("data_processed", showWarnings = FALSE)
saveRDS(dat, "data_processed/prad_master.rds")