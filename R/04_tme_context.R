# 04_tme_context.R
# Connect top module to tumour microenvironment

library(WGCNA)
library(dplyr)
library(tibble)
library(ggplot2)

# Load
tcga        <- readRDS("results/tcga.rds")
net         <- readRDS("results/wgcna_net.rds")
gex_filtered <- readRDS("results/gex_filtered.rds")
mae_tcga    <- curatedPCaData::getPCa("tcga")

# Module eigengenes
MEs <- net$MEs
MEs$sample_id <- rownames(MEs)

# CIBERSORT immune scores
cibersort <- mae_tcga[["cibersort"]] %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("sample_id")

# Merge
combined <- tcga$clin %>%
  inner_join(MEs, by = "sample_id") %>%
  inner_join(cibersort, by = "sample_id")

# Correlate MEgreen with immune cell types
immune_cols <- names(cibersort)[-1]

immune_cor <- sapply(immune_cols, function(cell) {
  ct <- cor.test(combined$MEgreen, combined[[cell]], 
                 method = "spearman", exact = FALSE)
  c(cor = ct$estimate, pval = ct$p.value)
})

immune_df <- as.data.frame(t(immune_cor)) %>%
  rownames_to_column("cell_type") %>%
  arrange(pval) %>%
  filter(pval < 0.05)

message("MEgreen correlates with:")
print(immune_df)

# ── Key figure: MEgreen by Gleason and immune context ────────────────────────
# Stratify by CD8 T cell high/low
cd8_col <- names(cibersort)[grepl("CD8", names(cibersort))][1]

if (!is.na(cd8_col)) {
  combined$cd8_group <- ifelse(
    combined[[cd8_col]] > median(combined[[cd8_col]], na.rm = TRUE),
    "CD8 high", "CD8 low"
  )
  
  p <- ggplot(combined, aes(x = gleason_binary, y = MEgreen, 
                            fill = cd8_group)) +
    geom_boxplot(outlier.shape = 21, alpha = 0.8) +
    scale_fill_manual(values = c("CD8 high" = "#2166ac", 
                                 "CD8 low"  = "#d6604d")) +
    labs(
      title    = "MEgreen eigengene by Gleason grade and CD8 T cell infiltration",
      subtitle = "TCGA prostate cancer cohort (n=229)",
      x        = "Gleason grade",
      y        = "MEgreen eigengene",
      fill     = "TME context"
    ) +
    theme_classic(base_size = 14)
  
  ggsave("figures/MEgreen_gleason_cd8.png", p, 
         width = 7, height = 5, dpi = 150)
  message("Figure saved.")
}

saveRDS(immune_df, "results/immune_cor.rds")
message("TME analysis complete.")