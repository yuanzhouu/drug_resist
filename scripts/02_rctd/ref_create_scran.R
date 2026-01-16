library(Seurat)
library(spacexr)
library(biomaRt)
library(stringr)
save_path <- "/data/xu_lab_projectsx/yuanzhou/drug_resist/results/02_rctd/ref"
dir.create(save_path, showWarnings = FALSE, recursive = TRUE)

# Load data
dat <- readRDS("/data/xu_lab_projectsx/Hasan/VisiumHD/rds/references/lam_lung.rds")

# remove Low-Quality cells
dat <- subset(dat, celltype_090523 != "Low-Quality")

# Create reference object
counts <- GetAssayData(dat, assay = "SoupX", layer = "counts")
cell_types <- dat$celltype_090523
clean_cell_types <- str_replace_all(cell_types, "[ /\\-]", "_")
clean_cell_types <- as.factor(clean_cell_types)
names(clean_cell_types) <- colnames(counts)

# subset genes
blocked_markers <- readRDS("/data/xu_lab_projectsx/yuanzhou/drug_resist/results/02_rctd/scran_markers/blocked_marker_vec.RDS")
counts <- counts[blocked_markers, ]

reference <- Reference(
    counts,
    clean_cell_types
  )

# save reference object
saveRDS(reference, file = paste0(save_path, "/ref_lam_lung_scran.rds"))