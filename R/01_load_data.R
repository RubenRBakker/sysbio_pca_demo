# 01_load_data.R
library(curatedPCaData)
library(MultiAssayExperiment)
library(dplyr)
library(tibble)

# Load cohorts
mae_tcga   <- getPCa("tcga")
mae_taylor <- getPCa("taylor")

# Helper function
prepare_cohort <- function(mae, gex_assay, cohort_name) {
  
  gex <- mae[[gex_assay]]
  
  clin <- colData(mae) %>%
    as.data.frame() %>%
    rownames_to_column("sample_id") %>%
    select(sample_id, gleason_grade, AR_activity, 
           tumor_purity_pathology) %>%
    mutate(
      gleason_binary = case_when(
        gleason_grade <= 6 ~ "low",
        gleason_grade >= 8 ~ "high",
        TRUE ~ NA_character_
      ),
      gleason_binary = factor(gleason_binary, levels = c("low", "high")),
      cohort = cohort_name
    ) %>%
    filter(!is.na(gleason_binary))
  
  common_samples <- intersect(clin$sample_id, colnames(gex))
  clin <- clin %>% filter(sample_id %in% common_samples)
  gex  <- gex[, common_samples]
  
  message(cohort_name, ": ", ncol(gex), " samples, ",
          sum(clin$gleason_binary == "low"), " low, ",
          sum(clin$gleason_binary == "high"), " high Gleason")
  
  list(gex = gex, clin = clin)
}

# Prepare cohorts
tcga   <- prepare_cohort(mae_tcga,   "gex.rsem.log", "TCGA")
taylor <- prepare_cohort(mae_taylor, "gex.rma",      "Taylor")

# Save
saveRDS(tcga,   "results/tcga.rds")
saveRDS(taylor, "results/taylor.rds")

message("Done.")