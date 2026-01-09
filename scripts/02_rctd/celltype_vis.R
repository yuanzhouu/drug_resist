library(Seurat)
library(ggplot2)
library(patchwork)
# save path
save_path <- "/data/xu_lab_projectsx/yuanzhou/drug_resist/results/01_celltype_vis"
dir.create(save_path, recursive = TRUE, showWarnings = FALSE)

# load data
d4 <- readRDS("/data/xu_lab_projectsx/Hasan/VisiumHD/rds/all samples/em_d4hi_SCT_RNA.rds")
d18 <- readRDS("/data/xu_lab_projectsx/Hasan/VisiumHD/rds/all samples/em_d18hi_SCT_RNA.rds")
d32 <- readRDS("/data/xu_lab_projectsx/Hasan/VisiumHD/rds/all samples/em_d32hi_SCT_RNA.rds")  
d44 <- readRDS("/data/xu_lab_projectsx/Hasan/VisiumHD/rds/all samples/em_d44hi_SCT_RNA.rds")
# Define celltypes for comparison
cell_types <- c("LAMCORE-1", "LAMCORE-2", "LAMCORE-3")

# Define samples
samples <- list(
  "d4" = d4,
  "d18" = d18,
  "d32" = d32,
  "d44" = d44
)

# Visualize 4 samples x 3 celltypes = 12 panels
# Each row represents one sample, each column represents one celltype
# Each celltype shown in red, others in gray
plot_list <- list()

for (sample_idx in seq_along(samples)) {
  sample_name <- names(samples)[sample_idx]
  sample_data <- samples[[sample_idx]]
  
  for (celltype_idx in seq_along(cell_types)) {
    celltype <- cell_types[celltype_idx]
    
    # Create binary column: target celltype vs others
    sample_data$highlight_celltype <- ifelse(sample_data$first_type == celltype, celltype, "Other")
    # Set factor levels to ensure correct color assignment
    sample_data$highlight_celltype <- factor(sample_data$highlight_celltype, levels = c("Other", celltype))
    
    # Create plot with colors - vector order must match factor level order
    # levels: "Other" (position 1) = gray, celltype (position 2) = red
    # image.alpha = 0 removes the H&E image background
    p <- SpatialDimPlot(sample_data, 
                        group.by = "highlight_celltype",
                        image.alpha = 0,
                        pt.size.factor = 1.7,
                        shape = 22) +
      ggtitle(paste0(sample_name, " - ", celltype))
    
    # Force override by adding manual scales (will override any default scales)
    p <- p + 
      scale_fill_manual(values = c("gray", "red"), drop = FALSE)
    
    # Store plot in list (row-major order: sample1_celltype1, sample1_celltype2, sample1_celltype3, sample2_celltype1, ...)
    plot_list[[(sample_idx - 1) * length(cell_types) + celltype_idx]] <- p
  }
}

# Combine all plots into one image: 4 rows (samples) x 3 columns (celltypes)
combined_plot <- wrap_plots(plot_list, ncol = 3, nrow = 4)

# Save combined plot
ggsave(paste0(save_path, "/LAMCORE_123_comparison_4x3panel.png"), 
       plot = combined_plot, width = 15, height = 20)
