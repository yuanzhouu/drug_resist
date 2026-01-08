library(Seurat)
library(ggplot2)
library(patchwork)
save_path <- "/data/xu_lab_projectsx/yuanzhou/drug_resist/results/01_celltype_vis"

# load original spaceranger output
em <- Load10X_Spatial(data.dir = "/data/xu_lab_projectsx/Hasan/VisiumHD/4h_em/outs", bin.size = c(8))
em # 254420 samples, 18085 genes

# show image
SpatialFeaturePlot(em, features = "nCount_Spatial.008um") + theme(legend.position = "right")
ggsave(paste0(save_path, "/spatial_image_em.png"), width = 5, height = 3)

# count how many spots have nCount_Spatial.008um > 10
sum(em@meta.data$nCount_Spatial.008um >= 10) # 240454

# summary
summary(em@meta.data$nCount_Spatial.008um)

# read existing rds file
em_rds <- readRDS("/data/xu_lab_projectsx/Hasan/VisiumHD/rds/all samples/em_d4hi_SCT_RNA.rds")
em_rds # 253950 samples, 35731 genes, 480 down from 254420 
# count the number of spots with nCount_Spatial.008um = 0
sum(em@meta.data$nCount_Spatial.008um == 0) # 480
sum(em_rds@meta.data$nCount_Spatial.008um >= 10) # 240454
# read existing RCTD file
rctd_rds <- readRDS("/data/xu_lab_projectsx/Hasan/VisiumHD/rds/all_RCTD/RCTD_d4hi_SCT.rds")
rctd_rds44 <- readRDS("/data/xu_lab_projectsx/Hasan/VisiumHD/rds/all_RCTD/RCTD_d44.rds")
str(rctd_rds@spatialRNA,max.level=2) # 250380 samples, 3560 down from 253950, only 3053 of nUMI < 10

# compare nUMI, confirmed is SCT assay is used for RCTD
counts <- GetAssayData(em_rds, assay = "RNA", layer = "counts")
sct_counts <- GetAssayData(em_rds, assay = "SCT", layer = "counts")

i = 11
barcode <- colnames(em_rds)[i]
colSums(counts)[i]
colSums(sct_counts)[i]
rctd_rds@originalSpatialRNA@nUMI[barcode]

# Calculate mean-variance relationship for RNA and SCT assays
# Calculate mean and variance for each gene in RNA counts
gene_means_rna <- Matrix::rowMeans(counts)
gene_vars_rna <- apply(counts, 1, var)

# Calculate mean and variance for each gene in SCT counts
gene_means_sct <- Matrix::rowMeans(sct_counts)
gene_vars_sct <- apply(sct_counts, 1, var)

# Create data frames for plotting
df_rna <- data.frame(mean = gene_means_rna, variance = gene_vars_rna)
df_sct <- data.frame(mean = gene_means_sct, variance = gene_vars_sct)

# Create scatter plots
p1 <- ggplot(df_rna, aes(x = mean, y = variance)) +
  geom_point(alpha = 0.3, size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "RNA Assay", x = "Mean Expression (log10)", y = "Variance (log10)") +
  theme_bw()

p2 <- ggplot(df_sct, aes(x = mean, y = variance)) +
  geom_point(alpha = 0.3, size = 0.5) +
  scale_x_log10() +
  scale_y_log10() +
  labs(title = "SCT Assay", x = "Mean Expression (log10)", y = "Variance (log10)") +
  theme_bw()

# Combine plots
combined_plot <- p1 + p2 + plot_annotation(title = "Mean-Variance Relationship by Assay")

# Save plot
ggsave(paste0(save_path, "/mean_variance_comparison.png"), combined_plot, width = 10, height = 5, dpi = 300)

#


# why reference only has 200 samples
str(rctd_rds@reference,max.level=2)
str(rctd_rds44@reference,max.level=2)
dim(rctd_rds@reference@counts)
ref <- readRDS("/data/xu_lab_projectsx/Hasan/VisiumHD/rds/references/lamref.rds")
str(ref,max.level=2)

# test RCTD code

### Load in/preprocess your data, this might vary based on your file type
refdir <- system.file("extdata",'Reference/Vignette',package = 'spacexr') # directory for the reference
counts <- read.csv(file.path(refdir,"dge.csv")) # load in counts matrix
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
meta_data <- read.csv(file.path(refdir,"meta_data.csv")) # load in meta_data (barcodes, clusters, and nUMI)
cell_types <- meta_data$cluster; names(cell_types) <- meta_data$barcode # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- meta_data$nUMI; names(nUMI) <- meta_data$barcode # create nUMI named list

### Create the Reference object
reference <- Reference(counts, cell_types, nUMI)

datadir <- system.file("extdata",'SpatialRNA/Vignette',package = 'spacexr') # directory for sample Slide-seq dataset
counts <- read.csv(file.path(datadir,"MappedDGEForR.csv")) # load in counts matrix
coords <- read.csv(file.path(datadir,"BeadLocationsForR.csv"))
rownames(counts) <- counts[,1]; counts[,1] <- NULL # Move first column to rownames
rownames(coords) <- coords$barcodes; coords$barcodes <- NULL # Move barcodes to rownames
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI

### Create SpatialRNA object
puck <- SpatialRNA(coords, counts, nUMI)

myRCTD <- create.RCTD(puck, reference, max_cores = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
results <- myRCTD@results

str(myRCTD@originalSpatialRNA,max.level=2)
str(myRCTD@spatialRNA,max.level=2)
dim(myRCTD@reference@counts)

max(table(reference@cell_types)) + 1

coerce_deglam_reference <- function(old_reference) {
  return(Reference(old_reference@counts, old_reference@cell_types,
                   old_reference@nUMI, n_max_cells = max(table(old_reference@cell_types)) + 1,
                   min_UMI = 1))
}

deglam_reference <- coerce_deglam_reference(reference)

myRCTD1 <- create.RCTD(puck, reference, max_cores = 1)
dim(myRCTD1@reference@counts)