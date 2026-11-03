# 02_wgcna.R
# WGCNA co-expression network on TCGA discovery cohort

library(WGCNA)
library(dplyr)

options(stringsAsFactors = FALSE)
allowWGCNAThreads(nThreads = 14)

# Load data
tcga <- readRDS("results/tcga.rds")

# Transpose: samples x genes
gex <- t(tcga$gex)

message("Full dimensions: ", nrow(gex), " samples x ", ncol(gex), " genes")

# Keep top 5000 most variable genes
gene_var <- apply(gex, 2, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:5000]
gex_filtered <- gex[, top_genes]

message("After filtering: ", nrow(gex_filtered), " x ", ncol(gex_filtered))

# Check for bad samples and genes
gsg <- goodSamplesGenes(gex_filtered, verbose = 3)
if (!gsg$allOK) {
  gex_filtered <- gex_filtered[gsg$goodSamples, gsg$goodGenes]
  message("Removed bad samples/genes")
}

# Soft threshold selection
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(
  gex_filtered,
  powerVector = powers,
  verbose = 5,
  networkType = "signed"
)
saveRDS(sft, "results/sft.rds")
message("Recommended soft threshold: ", sft$powerEstimate)

# Set power
soft_power <- sft$powerEstimate
if (is.na(soft_power)) soft_power <- 12

# Build network
net <- blockwiseModules(
  gex_filtered,
  power             = soft_power,
  networkType       = "signed",
  TOMType           = "signed",
  minModuleSize     = 30,
  mergeCutHeight    = 0.25,
  numericLabels     = FALSE,
  pamRespectsDendro = FALSE,
  saveTOMs          = FALSE,
  verbose           = 3,
  nThreads          = 14
)

saveRDS(net,          "results/wgcna_net.rds")
saveRDS(gex_filtered, "results/gex_filtered.rds")

message("WGCNA complete. Modules: ", length(unique(net$colors)) - 1)