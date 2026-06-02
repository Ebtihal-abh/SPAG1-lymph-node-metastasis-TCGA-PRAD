# ============================================
# 17_GSEA_visualisation.R  
# Purpose:
#   Generate publication-ready GSEA figures:
#   (A) Dotplot — SPAG1 Q4 vs Q1 (Script 08)
#   (B) Dotplot — Clinically adjusted N1 vs N0 (Script 08b)
#   (C) Concordance bar chart — 15 shared pathways
#   Combined 3-panel manuscript figure
#
# Inputs (all pre-computed, no new analysis):
#   tables/08/FGSEA_Hallmark_SPAG1_Q4_vs_Q1_multilevel.csv
#   tables/08b/FGSEA_Hallmark_N1_vs_N0_adjTGleason_TSS.csv
#   tables/08b/Pathway_overlap_N1_vs_SPAG1.csv
#
# Outputs:
#   figures/17/Fig_GSEA_SPAG1_dotplot.pdf
#   figures/17/Fig_GSEA_N1_dotplot.pdf
#   figures/17/Fig_GSEA_concordance.pdf
#   figures/17/Fig_GSEA_combined_panel.pdf
#   figures/17/Fig_GSEA_combined_panel.png
#
# ============================================

library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)


# ── 1. Output directory ───────────────────────────────────────────────────────
dir.create("figures/17", recursive = TRUE, showWarnings = FALSE)

# ── 2. Input file paths ───────────────────────────────────────────────────────
spag1_file <- "tables/08/FGSEA_Hallmark_SPAG1_Q4_vs_Q1_multilevel.csv"
n1_file      <- "tables/08b/FGSEA_Hallmark_N1_vs_N0_adjTGleason_TSS.csv"
overlap_file <- "tables/08b/Pathway_overlap_N1_vs_SPAG1.csv"

for (f in c(spag1_file, n1_file, overlap_file)) {
  if (!file.exists(f)) {
    stop("Required input file not found: ", f,
         "\nRun Scripts 08 and 08b first.")
  }
}

# ── 3. Load data ──────────────────────────────────────────────────────────────
spag1_gsea <- read.csv(spag1_file)
n1_gsea    <- read.csv(n1_file)
overlap    <- read.csv(overlap_file)

cat("SPAG1 GSEA loaded:", nrow(spag1_gsea), "pathways\n")
cat("N1 GSEA loaded:   ", nrow(n1_gsea),    "pathways\n")
cat("Overlap table:    ", nrow(overlap),     "pathways\n")

# Verify required columns
stopifnot(all(c("pathway", "NES", "padj") %in% colnames(spag1_gsea)))
stopifnot(all(c("pathway", "NES", "padj") %in% colnames(n1_gsea)))
stopifnot(all(c("pathway", "NES_N1", "NES_SPAG1", "Concordant")
              %in% colnames(overlap)))

# ── 4. Clean pathway names ────────────────────────────────────────────────────
clean_name <- function(x) {
  x <- sub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(tolower(x))
  
  # Restore acronyms broken by toTitleCase
  x <- gsub("\\bDna\\b",    "DNA",    x)
  x <- gsub("\\bUv\\b",     "UV",     x)
  x <- gsub("\\bMtorc1\\b", "mTORC1", x)
  x <- gsub("\\bNfkb\\b",   "NFkB",   x)
  x <- gsub("\\bIl6\\b",    "IL6",    x)
  x <- gsub("\\bJak\\b",    "JAK",    x)
  x <- gsub("\\bStat3\\b",  "STAT3",  x)
  x <- gsub("\\bStat\\b",   "STAT",   x)
  x <- gsub("\\bKras\\b",   "KRAS",   x)
  x <- gsub("\\bE2f\\b",    "E2F",    x)
  x <- gsub("\\bG2m\\b",    "G2M",    x)
  x <- gsub("\\bTgf\\b",    "TGF",    x)
  x <- gsub("\\bTnfa\\b",   "TNFa",   x)  
  x
}

spag1_gsea$pathway_clean <- clean_name(spag1_gsea$pathway)
n1_gsea$pathway_clean    <- clean_name(n1_gsea$pathway)
overlap$pathway_clean    <- clean_name(overlap$pathway)

# ── 5. Select top pathways for dotplots ───────────────────────────────────────
# Top 10 enriched (NES > 0) and top 10 depleted (NES < 0)
# from significant pathways only (padj < 0.05), ordered by padj

select_top <- function(df, n_each = 10) {
  enriched <- df %>%
    filter(padj < 0.05, NES > 0) %>%
    arrange(padj) %>%
    slice_head(n = n_each)
  
  depleted <- df %>%
    filter(padj < 0.05, NES < 0) %>%
    arrange(padj) %>%
    slice_head(n = n_each)
  
  bind_rows(enriched, depleted) %>%
    arrange(NES)
}

spag1_plot <- select_top(spag1_gsea, n_each = 10)
n1_plot    <- select_top(n1_gsea,    n_each = 10)

cat("\nSPAG1 pathways selected for dotplot:", nrow(spag1_plot), "\n")
cat("N1 pathways selected for dotplot:   ", nrow(n1_plot), "\n")

# ── 6. Compute shared scales ──────────────────────────────────────────────────
# Both dotplots must use identical x limits, colour limits, and size limits
# so panels A and B are directly visually comparable

# Shared NES limit — drives both x-axis and colour scale
nes_limit <- max(abs(c(spag1_plot$NES, n1_plot$NES)), na.rm = TRUE)
nes_limit <- ceiling(nes_limit * 10) / 10
cat("Shared NES limit (x-axis and colour):", nes_limit, "\n")

# Shared significance limit — drives dot size scale
sig_limit <- max(
  c(-log10(spag1_plot$padj), -log10(n1_plot$padj)),
  na.rm = TRUE
)
sig_limit <- ceiling(sig_limit)
cat("Shared -log10(padj) limit (dot size):", sig_limit, "\n")

# ── 7. Dotplot function ───────────────────────────────────────────────────────
# Produces a lollipop dotplot with:
#   - x position = NES
#   - dot colour = NES (blue-white-red)
#   - dot size = -log10(padj)
#   - shared scales passed as arguments so panels A and B are comparable
#   - NES colour legend above size legend (guides order)

make_dotplot <- function(df, title, subtitle, nes_lim, sig_lim) {
  
  # Order pathways by NES: most depleted at top, most enriched at bottom
  pathway_levels <- df$pathway_clean[order(df$NES)]
  df <- df %>%
    mutate(pathway_clean = factor(pathway_clean, levels = pathway_levels))
  
  ggplot(df, aes(x = NES, y = pathway_clean)) +
    
    # Lollipop stem from zero to dot
    geom_segment(
      aes(x = 0, xend = NES,
          y = pathway_clean, yend = pathway_clean),
      colour    = "grey80",
      linewidth = 0.4
    ) +
    
    # Dot: colour encodes NES, size encodes significance
    geom_point(aes(colour = NES, size = -log10(padj))) +
    
    # Colour scale: blue (depleted) to red (enriched), shared limits
    scale_colour_gradient2(
      low      = "#4575B4",
      mid      = "white",
      high     = "#D73027",
      midpoint = 0,
      limits   = c(-nes_lim, nes_lim),
      name     = "NES"
    ) +
    
    # Size scale: shared limits so dot sizes mean the same in both panels
    scale_size_continuous(
      range  = c(2, 7),
      limits = c(0, sig_lim),
      name   = expression(-log[10](padj))
    ) +
    
    # X-axis: shared limits
    scale_x_continuous(limits = c(-nes_lim, nes_lim)) +
    
    # Reference line at NES = 0
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey50") +
    
    labs(
      title    = title,
      subtitle = subtitle,
      x        = "Normalised Enrichment Score (NES)",
      y        = NULL
    ) +
    
    theme_classic(base_size = 11) +
    theme(
      plot.title         = element_text(face = "bold", size = 11),
      plot.subtitle      = element_text(size = 9, colour = "grey40"),
      axis.text.y        = element_text(size = 9),
      axis.text.x        = element_text(size = 9),
      legend.position    = "right",
      panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3)
    ) +
    
    # NES colour legend first, size legend second
    guides(
      colour = guide_colourbar(order = 1),
      size   = guide_legend(order = 2)
    )
}

# ── 8. Build dotplots ─────────────────────────────────────────────────────────
p_spag1 <- make_dotplot(
  df       = spag1_plot,
  title    = "SPAG1 High vs Low (Q4 vs Q1)",
  subtitle = "Hallmark GSEA | FDR < 0.05",
  nes_lim  = nes_limit,
  sig_lim  = sig_limit
)

p_n1 <- make_dotplot(
  df       = n1_plot,
  title    = "N1 vs N0",
  subtitle = "Adjusted for T stage, Gleason, TSS | Hallmark GSEA | FDR<0.05",
  nes_lim  = nes_limit,
  sig_lim  = sig_limit
)

# Save individual dotplots
ggsave(
  "figures/17/Fig_GSEA_SPAG1_dotplot.pdf",
  p_spag1,
  width  = 8,
  height = 6,
  device = "pdf"
)

ggsave(
  "figures/17/Fig_GSEA_N1_dotplot.pdf",
  p_n1,
  width  = 8,
  height = 6,
  device = "pdf"
)

cat("Individual dotplots saved\n")

# ── 9. Concordance bar chart ──────────────────────────────────────────────────
# Shows only the concordant overlapping pathways
# Paired horizontal bars: NES from SPAG1 analysis vs NES from N1 analysis
# Both analyses shown side by side for direct comparison

concordant <- overlap %>%
  filter(Concordant == TRUE)

cat("Concordant pathways for bar chart:", nrow(concordant), "\n")

# Reshape to long format for grouped bar chart
concordant_long <- concordant %>%
  select(pathway_clean, NES_N1, NES_SPAG1) %>%
  pivot_longer(
    cols      = c(NES_N1, NES_SPAG1),
    names_to  = "Analysis",
    values_to = "NES"
  ) %>%
  mutate(
    Analysis = case_match(Analysis,
                          "NES_N1"    ~ "N1 vs N0 (adj.)",
                          "NES_SPAG1" ~ "SPAG1 High vs Low"
    )
  )

pathway_order_conc <- concordant %>%
  arrange(NES_N1) %>%
  pull(pathway_clean)

concordant_long <- concordant_long %>%
  mutate(pathway_clean = factor(pathway_clean, levels = pathway_order_conc))

p_concordance <- ggplot(
  concordant_long,
  aes(x = NES, y = pathway_clean, fill = Analysis)
) +
  geom_col(
    position = position_dodge(width = 0.7),
    width    = 0.6,
    alpha    = 0.85
  ) +
  scale_fill_manual(
    values = c(
      "N1 vs N0 (adj.)"                 = "#2C7BB6",
      "SPAG1 High vs Low" = "#D7191C"
    )
  ) +
  geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
  labs(
    title = "Concordant transcriptional programmes in SPAG1-high tumours and node-positive disease",
    subtitle = paste0(
      nrow(concordant),
      " pathways concordant in direction (all FDR<0.05 in both analyses)"
    ),
    x    = "Normalised Enrichment Score (NES)",
    y    = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title         = element_text(face = "bold", size = 11),
    plot.subtitle      = element_text(size = 9, colour = "grey40"),
    axis.text.y        = element_text(size = 9),
    legend.position    = "bottom",
    legend.text        = element_text(size = 9),
    panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3)
  )

ggsave(
  "figures/17/Fig_GSEA_concordance.pdf",
  p_concordance,
  width  = 9,
  height = 6,
  device = "pdf"
)

cat("Concordance plot saved\n")

# ── 10. Combined 3-panel manuscript figure ────────────────────────────────────
# Layout: panels A and B side by side (top row)
#         panel C full width (bottom row)
# Tag levels A, B, C added 

combined <- (p_spag1 | p_n1) / p_concordance +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 12)
    )
  ) +
  plot_layout(heights = c(1, 0.85))

ggsave(
  "figures/17/Fig_GSEA_combined_panel.pdf",
  combined,
  width  = 16,
  height = 13,
  device = "pdf"
)

ggsave(
  "figures/17/Fig_GSEA_combined_panel.png",
  combined,
  width  = 16,
  height = 13,
  dpi    = 300
)

cat("Combined panel saved\n")

# ── 11. Final summary ─────────────────────────────────────────────────────────
cat("\n════ Script 17 Complete ════\n")
cat("Shared NES scale:           ±", nes_limit, "\n")
cat("Shared -log10(padj) max:    ",  sig_limit, "\n")
cat("Concordant pathways:        ",  nrow(concordant), "\n")
cat("\nOutputs saved to: figures/17/\n")
cat("  Fig_GSEA_SPAG1_dotplot.pdf     — Panel A\n")
cat("  Fig_GSEA_N1_dotplot.pdf        — Panel B\n")
cat("  Fig_GSEA_concordance.pdf       — Panel C\n")
cat("  Fig_GSEA_combined_panel.pdf    — Manuscript figure\n")
cat("  Fig_GSEA_combined_panel.png    — Preview\n")

