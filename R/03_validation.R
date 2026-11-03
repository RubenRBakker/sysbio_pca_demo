# 03_validation.R
# Module-Gleason correlation and Taylor validation

library(WGCNA)
library(dplyr)
library(ggplot2)
library(survival)
library(tibble)

# Load results
tcga        <- readRDS("results/tcga.rds")
taylor      <- readRDS("results/taylor.rds")
net         <- readRDS("results/wgcna_net.rds")
gex_filtered <- readRDS("results/gex_filtered.rds")

# ── Module eigengenes ─────────────────────────────────────────────────────────
MEs <- net$MEs
MEs$sample_id <- rownames(MEs)

# Merge with clinical
clin_tcga <- tcga$clin %>%
  inner_join(MEs, by = "sample_id")

# ── Correlate modules with Gleason ────────────────────────────────────────────
gleason_num <- as.numeric(clin_tcga$gleason_binary) - 1  # 0=low, 1=high

module_cols <- names(MEs)[grepl("^ME", names(MEs))]

cor_results <- sapply(module_cols, function(m) {
  ct <- cor.test(clin_tcga[[m]], gleason_num, method = "spearman")
  c(cor = ct$estimate, pval = ct$p.value)
})

cor_df <- as.data.frame(t(cor_results)) %>%
  rownames_to_column("module") %>%
  arrange(pval) %>%
  mutate(sig = pval < 0.05)

print(cor_df)

# ── Top modules ───────────────────────────────────────────────────────────────
top_modules <- cor_df %>% filter(sig) %>% pull(module)
message("Significant modules: ", paste(top_modules, collapse = ", "))

# ── Validate in Taylor ────────────────────────────────────────────────────────
# Project Taylor GEX onto TCGA modules
taylor_gex <- t(taylor$gex)

# Keep only genes present in both
common_genes <- intersect(colnames(gex_filtered), colnames(taylor_gex))
message("Common genes: ", length(common_genes))

taylor_gex_common  <- taylor_gex[, common_genes]
tcga_gex_common    <- gex_filtered[, common_genes]

# Calculate module eigengenes in Taylor using TCGA module assignments
module_assignments <- net$colors
common_assignments <- module_assignments[common_genes]

MEs_taylor <- moduleEigengenes(taylor_gex_common, common_assignments)$eigengenes
MEs_taylor$sample_id <- rownames(MEs_taylor)

clin_taylor <- taylor$clin %>%
  inner_join(MEs_taylor, by = "sample_id")

gleason_taylor <- as.numeric(clin_taylor$gleason_binary) - 1

# Validate top modules
val_results <- sapply(top_modules, function(m) {
  me_col <- gsub("ME", "ME", m)
  if (!me_col %in% names(clin_taylor)) return(c(cor = NA, pval = NA))
  ct <- cor.test(clin_taylor[[me_col]], gleason_taylor, method = "spearman")
  c(cor = ct$estimate, pval = ct$p.value)
})

val_df <- as.data.frame(t(val_results)) %>%
  rownames_to_column("module")

message("Validation results:")
print(val_df)

# Save
saveRDS(cor_df,  "results/module_gleason_cor.rds")
saveRDS(val_df,  "results/validation_results.rds")
saveRDS(MEs_taylor, "results/MEs_taylor.rds")

message("Validation complete.")