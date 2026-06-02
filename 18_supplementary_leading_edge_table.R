# ============================================
# 18_supplementary_leading_edge_table.R
# Purpose:
#   Generate Supplementary Table — leading edge genes
#   from the two primary GSEA analyses reported in the paper:
#   (1) SPAG1 Q4 vs Q1 (Script 08)
#   (2) N1 vs N0, clinically adjusted (Script 08b)
#
# Format: one row per gene per pathway (long format)
#   - Readable and sortable
#   - Includes NES, padj, direction for each pathway
#   - Significant pathways only (FDR < 0.05)
#
# Output:
#   tables/18/Supplementary_Table_Leading_Edge_Genes.csv
#   tables/18/Supplementary_Table_Leading_Edge_Genes_wide.csv
#   tables/18/Supplementary_Table_Concordant_Pathways_Leading_Edge.csv
# ============================================

library(dplyr)
library(tidyr)

dir.create("tables/18", recursive = TRUE, showWarnings = FALSE)

# ── 1. Load GSEA results ──────────────────────────────────────────────────────
spag1_adj <- read.csv(
  "tables/08/FGSEA_Hallmark_SPAG1_Q4_vs_Q1_multilevel.csv"
)
n1_adj <- read.csv(
  "tables/08b/FGSEA_Hallmark_N1_vs_N0_adjTGleason_TSS.csv"
)

cat("SPAG1 pathways loaded:", nrow(spag1_adj), "\n")
cat("N1 adj pathways loaded:   ", nrow(n1_adj), "\n")

# ── 2. Clean pathway names ────────────────────────────────────────────────────
clean_name <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  x <- gsub("\\bDna\\b",    "DNA",    x)
  x <- gsub("\\bUv\\b",     "UV",     x)
  x <- gsub("\\bMtorc1\\b", "mTORC1", x)
  x <- gsub("\\bNfkb\\b",   "NFkB",   x)
  x <- gsub("\\bIl6\\b",    "IL6",    x)
  x <- gsub("\\bJak\\b",    "JAK",    x)
  x <- gsub("\\bStat3\\b",  "STAT3",  x)
  x <- gsub("\\bKras\\b",   "KRAS",   x)
  x <- gsub("\\bE2f\\b",    "E2F",    x)
  x <- gsub("\\bG2m\\b",    "G2M",    x)
  x <- gsub("\\bTgf\\b",    "TGF",    x)
  x <- gsub("\\bStat\\b", "STAT", x)
  x <- gsub("\\bTnfa\\b",   "TNFa",   x)
  x
}

# ── 3. Function: expand leading edge to long format ───────────────────────────
# Takes a GSEA results dataframe and analysis label
# Returns one row per gene per significant pathway

expand_leading_edge <- function(df, analysis_label) {
  
  df %>%
    filter(padj < 0.05) %>%
    mutate(
      pathway_clean = clean_name(pathway),
      direction     = ifelse(NES > 0, "Enriched", "Depleted"),
      NES           = round(NES, 3),
      padj          = signif(padj, 3),
      analysis      = analysis_label
    ) %>%
    select(analysis, pathway_clean, direction, NES, padj, leadingEdge) %>%
    
    # Split semicolon-delimited leading edge into one row per gene
    mutate(gene = strsplit(leadingEdge, ";")) %>%
    unnest(gene) %>%
    mutate(gene = trimws(gene)) %>%
    filter(gene != "") %>%
    select(analysis, pathway_clean, direction, NES, padj, gene) %>%
    arrange(analysis, direction, padj, gene)
}

# ── 4. Expand both analyses ───────────────────────────────────────────────────
spag1_long <- expand_leading_edge(
  spag1_adj,
  "SPAG1 High vs Low"
)

n1_long <- expand_leading_edge(
  n1_adj,
  "N1 vs N0 (clinically adjusted)"
)

cat("\nSPAG1 long table rows:", nrow(spag1_long), "\n")
cat("N1 long table rows:   ", nrow(n1_long), "\n")

# ── 5. Combine into single supplementary table ───────────────────────────────
supp_table <- bind_rows(spag1_long, n1_long)

cat("Combined table rows:", nrow(supp_table), "\n")
cat("Columns:", paste(colnames(supp_table), collapse = ", "), "\n")

# ── 6. Rename columns for publication ────────────────────────────────────────
supp_table <- supp_table %>%
  rename(
    Analysis  = analysis,
    Pathway   = pathway_clean,
    Direction = direction,
    FDR       = padj,
    Gene      = gene
  )

# ── 7. Save long-format table ─────────────────────────────────────────────────
write.csv(
  supp_table,
  "tables/18/Supplementary_Table_Leading_Edge_Genes.csv",
  row.names = FALSE
)
cat("\nLong-format table saved\n")

# ── 8. Pathway summary table ──────────────────────────────────────────────────
# One row per pathway useful as a quick-reference header table
# to accompany the full leading edge table

pathway_summary <- supp_table %>%
  group_by(Analysis, Pathway, Direction, NES, FDR) %>%
  summarise(
    n_leading_edge_genes = n(),
    leading_edge_genes   = paste(sort(Gene), collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(Analysis, Direction, FDR)

write.csv(
  pathway_summary,
  "tables/18/Supplementary_Table_Pathway_Summary.csv",
  row.names = FALSE
)
cat("Pathway summary table saved\n")

# ── 9. Concordance subset ─────────────────────────────────────────────────────
# Leading edge genes specifically for the 15 concordant pathways
# This is the most interpretively important subset for the paper

overlap_file <- "tables/08b/Pathway_overlap_N1_vs_SPAG1.csv"
if (file.exists(overlap_file)) {
  
  overlap <- read.csv(overlap_file)
  concordant_pathways <- overlap %>%
    filter(Concordant == TRUE) %>%
    pull(pathway)
  
  concordant_pathways_clean <- clean_name(concordant_pathways)
  
  concordant_genes <- supp_table %>%
    filter(Pathway %in% concordant_pathways_clean)
  
  cat("\nConcordant pathway leading edge rows:", nrow(concordant_genes), "\n")
  cat("Concordant pathways covered:         ",
      length(unique(concordant_genes$Pathway)), "\n")
  
  write.csv(
    concordant_genes,
    "tables/18/Supplementary_Table_Concordant_Pathways_Leading_Edge.csv",
    row.names = FALSE
  )
  cat("Concordant pathway table saved\n")
  
} else {
  cat("Overlap file not found — skipping concordant subset\n")
}

# ── 10. Console preview ───────────────────────────────────────────────────────
cat("\n── Pathway summary preview ──\n")
print(as.data.frame(
  pathway_summary %>%
    select(Analysis, Pathway, Direction, NES, FDR, n_leading_edge_genes)
))

# ── 11. Summary ───────────────────────────────────────────────────────────────
cat("\n════ Script 18 Complete ════\n")
cat("Outputs saved to: tables/18/\n\n")
cat("  Supplementary_Table_Leading_Edge_Genes.csv\n")
cat("  — One row per gene per pathway\n")
cat("  — Both analyses combined\n")
cat("  — Total rows:", nrow(supp_table), "\n\n")
cat("  Supplementary_Table_Pathway_Summary.csv\n")
cat("  — One row per pathway with all leading edge genes in one cell\n")
cat("  — Total pathways:", nrow(pathway_summary), "\n\n")
cat("  Supplementary_Table_Concordant_Pathways_Leading_Edge.csv\n")
cat("  — Leading edge genes for the", 
    length(unique(concordant_genes$Pathway)), 
    "concordant pathways only\n")
cat("  — Highest interpretive priority for the paper\n")