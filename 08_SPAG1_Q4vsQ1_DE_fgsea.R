# 08_SPAG1_Q4vsQ1_DE_fgsea.R
# ============================================================
# Study:   TCGA-PRAD RNA-seq cohort (Illumina HiSeqV2)
# Purpose: Identify transcriptional programmes associated with
#          high vs low SPAG1 expression (Q4 vs Q1)
#          
#
# Inputs:
#   data_processed/prad_master_features.rds
#   data_processed/prad_expression_matrix.rds
#
# Outputs:
#   tables/08/DE_SPAG1_Q4_vs_Q1_limma.csv
#   tables/08/FGSEA_Hallmark_SPAG1_Q4_vs_Q1_multilevel.csv
# ============================================================

library(dplyr)
library(matrixStats)
library(limma)
library(msigdbr)
library(fgsea)
source("r/00_utils.R")

dir.create("tables/08", recursive = TRUE, showWarnings = FALSE)


# ---- Load data ----
dat <- readRDS("data_processed/prad_master_features.rds")
mat <- readRDS("data_processed/prad_expression_matrix.rds")

# ---- Alignment checks ----
stopifnot(nrow(dat) == ncol(mat))
stopifnot(all(dat$patient_id == substr(colnames(mat), 1, 12)))

# ---- Define Q1/Q4 groups on SPAG1 ----
q <- quantile(dat$SPAG1, probs = c(0.25, 0.75), na.rm = TRUE)
q1 <- unname(q[[1]])
q3 <- unname(q[[2]])

dat <- dat %>%
  mutate(
    SPAG1_group_q = case_when(
      SPAG1 <= q1 ~ "Low_Q1",
      SPAG1 >= q3 ~ "High_Q4",
      TRUE ~ NA_character_
    )
  )

keep <- !is.na(dat$SPAG1_group_q)
dat_sub <- dat[keep, , drop = FALSE]
mat_sub <- mat[, keep, drop = FALSE]
stopifnot(nrow(dat_sub) == ncol(mat_sub))

cat("Samples kept (Q1+Q4):", nrow(dat_sub), "\n")
print(table(dat_sub$SPAG1_group_q))

# ---- Common gene filter: remove zero-variance genes ----
gene_var <- matrixStats::rowVars(mat_sub)
keep_genes <- gene_var > 0
mat_sub_f <- mat_sub[keep_genes, , drop = FALSE]

cat("Genes before:", nrow(mat_sub), "\n")
cat("Genes after var>0 filter:", nrow(mat_sub_f), "\n")

# ---- Hallmark gene sets (prepare once) ----
hallmark <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)

# =========================================================
# A) UNADJUSTED limma: ~ group
# =========================================================
cat("\n=== A) Unadjusted limma: ~ group ===\n")

group <- factor(dat_sub$SPAG1_group_q, levels = c("Low_Q1", "High_Q4"))
design_unadj <- model.matrix(~ group)

fit_unadj <- lmFit(mat_sub_f, design_unadj)
fit_unadj <- eBayes(fit_unadj)

de_unadj <- topTable(fit_unadj, coef = "groupHigh_Q4", number = Inf, sort.by = "P")
write.csv(de_unadj, "tables/08/DE_SPAG1_Q4_vs_Q1_limma.csv", row.names = TRUE)
cat("Saved: tables/08/DE_SPAG1_Q4_vs_Q1_limma.csv\n")

# fgsea ranks from unadjusted model
ranks_unadj <- de_unadj$t
names(ranks_unadj) <- rownames(de_unadj)
set.seed(123)
ranks_unadj <- ranks_unadj + rnorm(length(ranks_unadj), 0, 1e-10)
ranks_unadj <- sort(ranks_unadj, decreasing = TRUE)
fg_unadj <- fgseaMultilevel(pathways = pathways, stats = ranks_unadj) %>%
  arrange(padj)

fg_unadj_out <- fg_unadj %>%
  mutate(leadingEdge = vapply(leadingEdge, function(x) paste(x, collapse = ";"), character(1)))

write.csv(fg_unadj_out,
          "tables/08/FGSEA_Hallmark_SPAG1_Q4_vs_Q1_multilevel.csv",
          row.names = FALSE)
cat("Saved: tables/08/FGSEA_Hallmark_SPAG1_Q4_vs_Q1_multilevel.csv\n")


# ---- Summary: top pathways ----
cat("\nTop significant pathways (unadjusted, FDR < 0.05):\n")
print(fg_unadj_out %>% filter(padj < 0.05) %>%
        select(pathway, NES, padj) %>%
        arrange(padj),
      row.names = FALSE)
cat("\nDONE.\n")