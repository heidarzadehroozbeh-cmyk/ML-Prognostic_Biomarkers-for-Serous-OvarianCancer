# scripts/02_download_geo.R
# Purpose: Download GEO GSEMatrix objects (mRNA) and save:
#   1) RDS of full GSEMatrix list
#   2) CSV of pData (metadata) for the first ExpressionSet (eset1)
#
# Note: If the RDS file already exists, the dataset will NOT be downloaded again.
#
# Run:
#   source("scripts/02_download_geo.R")

source("scripts/01_params.R")

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
})

message("Project root: ", getwd())
message("Saving GEO objects to: ", DIR_GEO)
message("Saving metadata to: ", DIR_META)

download_gse <- function(gse_id) {
  
  out_rds <- file.path(DIR_GEO, paste0(gse_id, "_GSEMatrix_list.rds"))
  
  # ---- Skip download if file already exists ----
  if (file.exists(out_rds)) {
    message("\n--- ", gse_id, " already exists. Skipping download. ---")
    return(NULL)
  }
  
  message("\n--- Downloading ", gse_id, " ---")
  
  # GEOquery returns a list of ExpressionSets
  gse_list <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE)
  
  # Save the full list
  saveRDS(gse_list, out_rds)
  
  # Write metadata of first ExpressionSet
  eset1 <- gse_list[[1]]
  pheno <- pData(eset1)
  
  out_csv <- file.path(DIR_META, paste0(gse_id, "_pData_eset1.csv"))
  write.csv(pheno, out_csv, row.names = TRUE)
  
  # Log
  message("Saved: ", out_rds)
  message("Saved: ", out_csv)
  message("Platforms in this GSE: ", paste(annotation(eset1), collapse = ", "))
  message("Samples in eset1: ", ncol(eset1))
}

# Download all mRNA datasets
for (id in GSE_mrna) {
  download_gse(id)
}

message("\nAll GEO downloads completed successfully.")
