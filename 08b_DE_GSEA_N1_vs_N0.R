# =============================================================================
# Script 08b: Differential Expression & GSEA — N0 vs N1
# =============================================================================
# Study:      TCGA-PRAD RNA-seq cohort (Illumina HiSeqV2)
# Comparison: N1 (node-positive) vs N0 (node-negative)
# Purpose:    Characterize the transcriptional programs associated with
#             lymph node metastasis, and assess overlap with SPAG1-associated
#             pathways from Script 08 to evaluate mechanistic coherence.
#
# Approach:
#   A) Unadjusted DE + GSEA: N1 vs N0 as sole predictor
#   B) Clinically adjusted DE + GSEA: N1 vs N0 adjusted for T stage and
#      Gleason score (removes confounding by known clinical predictors,
#      isolating the molecular signal specific to nodal spread)
#   C) Overlap analysis: compare top enriched pathways with Script 08
#      (SPAG1 Q4 vs Q1) to assess mechanistic coherence
#   
# NOTE on design:
#   N1 = 1 (positive), N0 = 0 (reference)
#   Positive NES in GSEA = pathway enriched in N1
#   Negative NES = pathway enriched in N0
#
# Inputs:
#   data_processed/prad_master_features.rds
#   data_processed/prad_expression_matrix.rds
#   tables/08/FGSEA_Hallmark_SPAG1_Q4_vs_Q1_multilevel.csv
#
# Outputs:
#   tables/08b/TSS_N1_rate_distribution.csv
#   tables/08b/TSS_collapsed_N1_balance.csv
#   tables/08b/DE_N1_vs_N0_unadjusted_limma.csv
#   tables/08b/FGSEA_Hallmark_N1_vs_N0_unadjusted.csv
#   tables/08b/DE_N1_vs_N0_adjTGleason_TSS_limma.csv
#   tables/08b/FGSEA_Hallmark_N1_vs_N0_adjTGleason_TSS.csv
#   tables/08b/Pathway_overlap_N1_vs_SPAG1.csv
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(matrixStats)
  library(limma)
  library(msigdbr)
  library(fgsea)
})

source("r/00_utils.R")

dir.create("tables/08b",  recursive = TRUE, showWarnings = FALSE)
dir.create("figures/08b", recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 0. LOAD DATA
# =============================================================================

dat <- readRDS("data_processed/prad_master_features.rds")
mat <- readRDS("data_processed/prad_expression_matrix.rds")

# Alignment check
stopifnot(nrow(dat) == ncol(mat))
stopifnot(all(dat$patient_id == substr(colnames(mat), 1, 12)))

# =============================================================================
# 1. DEFINE N STAGE GROUPS
# =============================================================================

dat <- dat %>%
  mutate(
    N_binary = case_when(
      tolower(as.character(N_group2)) == "n1" ~ 1L,
      tolower(as.character(N_group2)) == "n0" ~ 0L,
      TRUE ~ NA_integer_
    )
  )

# Restrict to patients with known N stage
keep_N <- !is.na(dat$N_binary)
dat_N  <- dat[keep_N, , drop = FALSE]
mat_N  <- mat[, keep_N, drop = FALSE]

stopifnot(nrow(dat_N) == ncol(mat_N))

cat("=============================================================\n")
cat(" SAMPLE SUMMARY — N stage analysis\n")
cat("=============================================================\n")
cat(sprintf(" Total with known N stage : %d\n", nrow(dat_N)))
cat(sprintf(" N0 (node-negative)       : %d\n", sum(dat_N$N_binary == 0)))
cat(sprintf(" N1 (node-positive)       : %d\n", sum(dat_N$N_binary == 1)))
cat(sprintf(" N1 rate                  : %.1f%%\n\n",
            100 * mean(dat_N$N_binary)))

# =============================================================================
# 2. GENE FILTERING — remove zero-variance genes
# =============================================================================

gene_var   <- matrixStats::rowVars(mat_N)
keep_genes <- gene_var > 0
mat_N_f    <- mat_N[keep_genes, , drop = FALSE]

cat(sprintf(" Genes before filter : %d\n", nrow(mat_N)))
cat(sprintf(" Genes after var > 0 : %d\n\n", nrow(mat_N_f)))

# =============================================================================
# 3. HALLMARK GENE SETS
# =============================================================================

hallmark  <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::select(gs_name, gene_symbol)
pathways  <- split(hallmark$gene_symbol, hallmark$gs_name)
cat(sprintf(" Hallmark pathways loaded: %d\n\n", length(pathways)))

# =============================================================================
# HELPER: run fgsea and format output
# =============================================================================

run_fgsea <- function(t_stats, gene_names, seed = 123) {
  ranks <- t_stats
  names(ranks) <- gene_names
  set.seed(seed)
  ranks <- ranks + rnorm(length(ranks), 0, 1e-10)  # break ties
  ranks <- sort(ranks, decreasing = TRUE)
  
  fg <- fgseaMultilevel(pathways = pathways, stats = ranks) %>%
    arrange(padj) %>%
    mutate(
      leadingEdge = vapply(leadingEdge,
                           function(x) paste(x, collapse = ";"),
                           character(1))
    )
  return(fg)
}

# =============================================================================
# BATCH EFFECT CHECK — is N stage confounded with tissue source site?
# =============================================================================

tss_n1 <- dat_N %>%
  group_by(tissue_source_site) %>%
  summarise(
    n_total = n(),
    n_N1    = sum(N_binary == 1, na.rm = TRUE),
    N1_rate = round(100 * mean(N_binary == 1, na.rm = TRUE), 1)
  ) %>%
  filter(n_total >= 5) %>%
  arrange(desc(N1_rate))

write.csv(tss_n1, "tables/08b/TSS_N1_rate_distribution.csv",
          row.names = FALSE)

tss_table  <- table(dat_N$tissue_source_site, dat_N$N_binary)
tss_table  <- tss_table[rowSums(tss_table) >= 5, ]
chi_result <- chisq.test(tss_table, simulate.p.value = TRUE, B = 2000)

cat(sprintf(" TSS batch check — N1 rate range: %.1f%% to %.1f%% (SD = %.1f%%)\n",
            min(tss_n1$N1_rate), max(tss_n1$N1_rate), sd(tss_n1$N1_rate)))
cat(sprintf(" Chi-square test (N stage vs TSS): p = %.4f\n",
            chi_result$p.value))
cat(sprintf(" Conclusion: TSS %s added as covariate in adjusted model.\n\n",
            ifelse(chi_result$p.value < 0.05, "WILL BE", "will NOT be")))
# =============================================================================
# A) UNADJUSTED DE + GSEA — N1 vs N0
# =============================================================================

cat("=============================================================\n")
cat(" A) Unadjusted: N1 vs N0\n")
cat("=============================================================\n")

N_factor      <- factor(dat_N$N_binary, levels = c(0, 1),
                        labels = c("N0", "N1"))
design_unadj  <- model.matrix(~ N_factor)

fit_unadj  <- lmFit(mat_N_f, design_unadj)
fit_unadj  <- eBayes(fit_unadj)

de_unadj <- topTable(fit_unadj, coef = "N_factorN1",
                     number = Inf, sort.by = "P")

write.csv(de_unadj,
          "tables/08b/DE_N1_vs_N0_unadjusted_limma.csv",
          row.names = TRUE)
cat(sprintf(" Saved: tables/08b/DE_N1_vs_N0_unadjusted_limma.csv\n"))
cat(sprintf(" Sig DE genes (FDR < 0.05): %d\n",
            sum(de_unadj$adj.P.Val < 0.05, na.rm = TRUE)))
cat(sprintf(" Top upregulated in N1: %s\n",
            paste(head(rownames(de_unadj[de_unadj$logFC > 0, ]), 5),
                  collapse = ", ")))
cat(sprintf(" Top downregulated in N1: %s\n\n",
            paste(head(rownames(de_unadj[de_unadj$logFC < 0, ]), 5),
                  collapse = ", ")))

fg_unadj <- run_fgsea(de_unadj$t, rownames(de_unadj))

write.csv(fg_unadj,
          "tables/08b/FGSEA_Hallmark_N1_vs_N0_unadjusted.csv",
          row.names = FALSE)
cat(" Saved: tables/08b/FGSEA_Hallmark_N1_vs_N0_unadjusted.csv\n")
cat(" Top enriched pathways in N1 (unadjusted):\n")
print(fg_unadj %>%
        filter(NES > 0, padj < 0.05) %>%
        select(pathway, NES, padj) %>%
        head(10),
      row.names = FALSE)
cat("\n")

# =============================================================================
# B) Clinically adjusted DE + GSEA: N1 vs N0 adjusted for T stage,
#      Gleason score, and tissue source site / TSS (removes confounding
#      by clinical predictors and center-level surgical practice variation,
#      isolating the molecular signal specific to nodal spread)
# =============================================================================

cat("=============================================================\n")
cat(" B) Adjusted: N1 vs N0 + T stage + Gleason + TSS\n")
cat("=============================================================\n")

model_vars_adj <- c("N_binary", "T_group2", "Gleason_group2")
keep_adj <- complete.cases(dat_N[, model_vars_adj])
cat(sprintf(" Dropped due to missing T/Gleason: %d\n", sum(!keep_adj)))

dat_adj <- dat_N[keep_adj, , drop = FALSE]
mat_adj <- mat_N_f[, keep_adj, drop = FALSE]

dat_adj$TSS_collapsed <- collapse_tss(dat_adj$tissue_source_site)

# Quick balance check — N1 distribution across TSS levels
balance_check <- dat_adj %>%
  group_by(TSS_collapsed) %>%
  summarise(
    n_total = n(),
    n_N1    = sum(N_binary == 1),
    N1_rate = round(100 * mean(N_binary == 1), 1)
  ) %>%
  arrange(desc(N1_rate))

write.csv(balance_check, 
          "tables/08b/TSS_collapsed_N1_balance.csv",
          row.names = FALSE)
print(balance_check)

cat(sprintf(" TSS levels after collapsing (>=25 samples): %d\n",
            nlevels(dat_adj$TSS_collapsed)))
cat(sprintf(" Adjusted-model sample size : %d\n", nrow(dat_adj)))
cat(sprintf(" N1 in adjusted sample      : %d\n\n",
            sum(dat_adj$N_binary == 1)))

dat_adj$T_group2       <- factor(dat_adj$T_group2)
dat_adj$Gleason_group2 <- factor(dat_adj$Gleason_group2)
dat_adj$N_factor       <- factor(dat_adj$N_binary,
                                 levels = c(0, 1),
                                 labels = c("N0", "N1"))

design_adj <- model.matrix(
  ~ TSS_collapsed + T_group2 + Gleason_group2 + N_factor,
  data = dat_adj
)

cat(sprintf(" Design matrix columns: %d\n\n", ncol(design_adj)))

fit_adj  <- lmFit(mat_adj, design_adj)
fit_adj  <- eBayes(fit_adj)

de_adj <- topTable(fit_adj, coef = "N_factorN1",
                   number = Inf, sort.by = "P")

write.csv(de_adj,
          "tables/08b/DE_N1_vs_N0_adjTGleason_TSS_limma.csv",
          row.names = TRUE)
cat(sprintf(" Saved: DE_N1_vs_N0_adjTGleason_TSS_limma.csv\n"))
cat(sprintf(" Sig DE genes (FDR < 0.05): %d\n",
            sum(de_adj$adj.P.Val < 0.05, na.rm = TRUE)))
cat(sprintf(" SPAG1 in results: logFC = %.3f, adj.P.Val = %.4f\n",
            de_adj["SPAG1", "logFC"],
            de_adj["SPAG1", "adj.P.Val"]))

fg_adj <- run_fgsea(de_adj$t, rownames(de_adj))

write.csv(fg_adj,
          "tables/08b/FGSEA_Hallmark_N1_vs_N0_adjTGleason_TSS.csv",
          row.names = FALSE)
cat(" Saved: FGSEA_Hallmark_N1_vs_N0_adjTGleason_TSS.csv\n")
cat(" Top enriched pathways in N1 (TSS-adjusted):\n")
print(fg_adj %>%
        filter(NES > 0, padj < 0.05) %>%
        select(pathway, NES, padj) %>%
        head(10),
      row.names = FALSE)
cat("\n")

# =============================================================================
# C) OVERLAP ANALYSIS — compare with Script 08 SPAG1 Q4 vs Q1 results
# =============================================================================

spag1_file <- "tables/08/FGSEA_Hallmark_SPAG1_Q4_vs_Q1_multilevel.csv"

cat("=============================================================\n")
cat(" C) Overlap with SPAG1 Q4 vs Q1 (Script 08)\n")
cat("=============================================================\n")

fg_adj_file <- "tables/08b/FGSEA_Hallmark_N1_vs_N0_adjTGleason_TSS.csv"
fg_adj <- read.csv(fg_adj_file) 

if (file.exists(spag1_file)) {
  fg_spag1 <- read.csv(spag1_file)
  
  # Significant pathways in each analysis (FDR < 0.05)
  sig_N1    <- fg_adj %>%
    filter(padj < 0.05) %>%
    pull(pathway)
  
  sig_SPAG1 <- fg_spag1 %>%
    filter(padj < 0.05) %>%
    pull(pathway)
  
  overlap <- intersect(sig_N1, sig_SPAG1)
  
  cat(sprintf(" Sig pathways in N1 analysis (FDR<0.05)    : %d\n",
              length(sig_N1)))
  cat(sprintf(" Sig pathways in SPAG1 analysis (FDR<0.05) : %d\n",
              length(sig_SPAG1)))
  cat(sprintf(" Overlapping pathways                       : %d\n\n",
              length(overlap)))
  
  if (length(overlap) > 0) {
    # Get NES direction from both analyses for overlapping pathways
    nes_N1 <- fg_adj %>%
      filter(pathway %in% overlap) %>%
      select(pathway, NES_N1 = NES, padj_N1 = padj)
    
    nes_SPAG1 <- fg_spag1 %>%
      filter(pathway %in% overlap) %>%
      select(pathway, NES_SPAG1 = NES, padj_SPAG1 = padj)
    
    overlap_tbl <- left_join(nes_N1, nes_SPAG1, by = "pathway") %>%
      mutate(
        Direction_N1    = ifelse(NES_N1 > 0, "Enriched in N1", "Enriched in N0"),
        Direction_SPAG1 = ifelse(NES_SPAG1 > 0, "Enriched in High-SPAG1",
                                 "Enriched in Low-SPAG1"),
        Concordant = (NES_N1 > 0 & NES_SPAG1 > 0) |
          (NES_N1 < 0 & NES_SPAG1 < 0)
      ) %>%
      arrange(desc(Concordant), padj_N1)
    
    write.csv(overlap_tbl,
              "tables/08b/Pathway_overlap_N1_vs_SPAG1.csv",
              row.names = FALSE)
    cat(" Saved: tables/08b/Pathway_overlap_N1_vs_SPAG1.csv\n\n")
    
    cat(" Overlapping pathways with direction:\n")
    print(overlap_tbl %>%
            select(pathway, NES_N1, NES_SPAG1, Concordant),
          row.names = FALSE)
    
    n_concordant <- sum(overlap_tbl$Concordant)
    cat(sprintf("\n Concordant direction (both point same way): %d / %d\n",
                n_concordant, length(overlap)))
    
    if (n_concordant == length(overlap)) {
      cat(" Interpretation: ALL overlapping pathways are concordant.\n")
      cat(" HIGH mechanistic coherence — SPAG1 biology aligns with N1 biology.\n")
    } else if (n_concordant / length(overlap) >= 0.7) {
      cat(" Interpretation: Most overlapping pathways are concordant.\n")
      cat(" GOOD mechanistic coherence.\n")
    } else {
      cat(" Interpretation: Mixed directions. Limited mechanistic coherence.\n")
    }
    
  } else {
    cat(" No overlapping significant pathways found.\n")
    cat(" SPAG1 and N1 transcriptional programs appear distinct.\n")
    cat(" Recommendation: do not add N0 vs N1 analysis to paper.\n")
  }
  
} else {
  cat(sprintf(" Script 08 results not found at: %s\n", spag1_file))
  cat(" Run Script 08 first, then rerun this overlap section.\n")
}

cat("\n=============================================================\n")
cat(" Script 08b complete.\n")
cat("=============================================================\n")