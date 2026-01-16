library(Seurat)
library(spacexr)

save_path <- "/data/xu_lab_projectsx/yuanzhou/drug_resist/results/02_rctd/query"
dir.create(save_path, showWarnings = FALSE, recursive = TRUE)

samples <- c("4h_em", "18h_em", "32h_em", "44h_em")

# Process each sample
for (sample in samples) {
  cat("Processing sample:", sample, "\n")
  
  # Load data
  dat <- readRDS(paste0("/data/xu_lab_projectsx/yuanzhou/drug_resist/results/01_qc/", sample, "/", sample, "_qc.rds"))
  
  # Create query object
  coords <- GetTissueCoordinates(dat)[,1:2]
  counts <- GetAssayData(dat, assay = "Spatial.008um", layer = "counts")
  rownames(counts) <- dat[["Spatial.008um"]][["ENSG"]][[1]]
  
  query <- SpatialRNA(coords, counts)
  
  # Save query object
  saveRDS(query, file = paste0(save_path, "/query_", sample, ".rds"))
  
  cat("Saved query object for", sample, "\n")
}
