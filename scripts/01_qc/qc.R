# This code is used to perform qc on the spatial transcriptomics data.

# qc removal criterion:
# 1. mitocondrial gene percentage > 25%
# 2. number of features 2.5% in the tail
# 3. number of counts 2.5% in the tail

save_path <- "/data/xu_lab_projectsx/yuanzhou/drug_resist/results/01_qc"

library(Seurat)
library(ggplot2)
# data loading

samples <- c("4h_em", "18h_em", "32h_em", "44h_em")

# Loop through all samples
for (sample in samples) {
  cat("\nProcessing sample:", sample, "\n")
  
  # create subfolder storing qc plots
  sample_path <- paste0("/data/xu_lab_projectsx/yuanzhou/drug_resist/results/01_qc/", sample)
  dir.create(sample_path, showWarnings = FALSE, recursive = TRUE)

  dat <- Load10X_Spatial(data.dir = paste0("/data/xu_lab_projectsx/Hasan/VisiumHD/", sample, "/outs"), bin.size = c(8))

  # add mitocondrial gene percentage to the data
  dat[["percent.mt"]] <- PercentageFeatureSet(dat, pattern = "^MT-")

  # Visualize QC metrics as a violin plot (before filtering)
  VlnPlot(dat, features = c("percent.mt", "nCount_Spatial.008um", "nFeature_Spatial.008um"), 
          ncol = 3, pt.size = 0)
  ggsave(paste0(sample_path, "/", sample, "_qc_violin_before.png"), width = 12, height = 4)

  # cell removal
  # Calculate thresholds for 2.5% tail (both lower and upper)
  nFeature_threshold_lower <- quantile(dat$nFeature_Spatial.008um, 0.025)
  nFeature_threshold_upper <- quantile(dat$nFeature_Spatial.008um, 0.975)
  nCount_threshold_lower <- quantile(dat$nCount_Spatial.008um, 0.025)
  nCount_threshold_upper <- quantile(dat$nCount_Spatial.008um, 0.975)

  # Filter cells based on QC criteria
  dat_filtered <- subset(dat, 
                         subset = percent.mt < 25 & 
                                  nFeature_Spatial.008um > nFeature_threshold_lower &
                                  nFeature_Spatial.008um < nFeature_threshold_upper &
                                  nCount_Spatial.008um > nCount_threshold_lower &
                                  nCount_Spatial.008um < nCount_threshold_upper)

  # Print filtering statistics and save to file
  log_file <- paste0(sample_path, "/", sample, "_qc_statistics.txt")
  sink(log_file, split = TRUE)
  cat("QC Filtering Statistics\n")
  cat("=======================\n\n")
  cat("Before filtering:", ncol(dat), "spots\n")
  cat("After filtering:", ncol(dat_filtered), "spots\n")
  cat("Removed:", ncol(dat) - ncol(dat_filtered), "spots\n")
  cat("\nFiltering thresholds:\n")
  cat("  Mitochondrial percentage threshold: < 25%\n")
  cat("  nFeature lower threshold (2.5%):", nFeature_threshold_lower, "\n")
  cat("  nFeature upper threshold (97.5%):", nFeature_threshold_upper, "\n")
  cat("  nCount lower threshold (2.5%):", nCount_threshold_lower, "\n")
  cat("  nCount upper threshold (97.5%):", nCount_threshold_upper, "\n")
  sink()

  # Print quantiles of nFeature_Spatial.008um
  nFeature_quantiles <- quantile(dat$nFeature_Spatial.008um, c(0.1, 0.3, 0.5, 0.7, 0.9))
  sink(log_file, append = TRUE, split = TRUE)
  cat("\nnFeature_Spatial.008um quantiles:\n")
  cat("  10%:", nFeature_quantiles[1], "\n")
  cat("  30%:", nFeature_quantiles[2], "\n")
  cat("  50%:", nFeature_quantiles[3], "\n")
  cat("  70%:", nFeature_quantiles[4], "\n")
  cat("  90%:", nFeature_quantiles[5], "\n")
  sink()

  # print quantiles of nCount_Spatial.008um
  nCount_quantiles <- quantile(dat$nCount_Spatial.008um, c(0.1, 0.3, 0.5, 0.7, 0.9))
  sink(log_file, append = TRUE, split = TRUE)
  cat("\nnCount_Spatial.008um quantiles:\n")
  cat("  10%:", nCount_quantiles[1], "\n")
  cat("  30%:", nCount_quantiles[2], "\n")
  cat("  50%:", nCount_quantiles[3], "\n")
  cat("  70%:", nCount_quantiles[4], "\n")
  cat("  90%:", nCount_quantiles[5], "\n")
  sink()
  
  # revisualization
  VlnPlot(dat_filtered, features = c("percent.mt", "nCount_Spatial.008um", "nFeature_Spatial.008um"), 
          ncol = 3, pt.size = 0)
  ggsave(paste0(sample_path, "/", sample, "_qc_violin_after.png"), width = 12, height = 4)

  # save RDS file, name it + _qc tag, save in sample_path
  saveRDS(dat_filtered, file = paste0(sample_path, "/", sample, "_qc.rds"))
  
  cat("Completed processing sample:", sample, "\n")
}