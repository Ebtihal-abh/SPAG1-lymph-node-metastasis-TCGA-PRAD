# Define and factor clinical grouping variables (T stage, N stage, Gleason)
# from TCGA-PRAD master dataset. Saves prad_master_grouped.rds.

library(dplyr)
library(stringr)

dat <- readRDS("data_processed/prad_master.rds")

# Clean T stage
dat <- dat %>%
  mutate(
    T_raw = toupper(as.character(pathologic_T)),
    T_raw = str_replace(T_raw, "^PT", "T"),
    T_group2 = case_when(
      str_detect(T_raw, "^T2") ~ "T2",
      str_detect(T_raw, "^T3|^T4") ~ "T3_4",
      TRUE ~ NA_character_
    )
  )

# Clean N stage
dat <- dat %>%
  mutate(
    N_raw = toupper(as.character(pathologic_N)),
    N_group2 = case_when(
      N_raw == "N0" ~ "N0",
      N_raw == "N1" ~ "N1",
      TRUE ~ NA_character_
    )
  )

dat <- dat %>%
  mutate(
    T_group2 = factor(T_group2, levels = c("T2", "T3_4")),
    N_group2 = factor(N_group2, levels = c("N0", "N1"))
  )


dat <- dat %>%
  mutate(
    gleason_num = as.numeric(gleason_score),
    Gleason_group2 = case_when(
      gleason_num >= 8 ~ "High",
      gleason_num <= 7 ~ "Low_intermediate",
      TRUE ~ NA_character_
    ),
    
    Gleason_group2 = factor(Gleason_group2,
                            levels = c("Low_intermediate", "High"))
  )

saveRDS(dat, "data_processed/prad_master_grouped.rds")


