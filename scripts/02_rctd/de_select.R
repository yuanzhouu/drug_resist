library(Seurat)
library(spacexr)
library(scran)
library(stringr)
library(ggplot2)

save_path <- '/data/xu_lab_projectsx/yuanzhou/drug_resist/results/02_rctd/scran_markers'
# check if RCTD object stores DE information

RCTD_4h_em <- readRDS("/data/xu_lab_projectsx/yuanzhou/drug_resist/results/02_rctd/test_results/4h_em/RCTD_4h_em.rds")
de_genes <- RCTD_4h_em@internal_vars$gene_list_reg
#platform_genes <- RCTD_4h_em@internal_vars$gene_list_bulk

# load single cell data and calculate DE genes
dat <- readRDS("/data/xu_lab_projectsx/Hasan/VisiumHD/rds/references/lam_lung.rds")

# remove Low-Quality cells
dat <- subset(dat, celltype_090523 != "Low-Quality")

# set assay to SoupX
DefaultAssay(dat) <- "SoupX"
# lognormalize SoupX assay
dat <- NormalizeData(dat)
# extract lognormalized counts
counts <- GetAssayData(dat, assay = "SoupX", layer = "data")
cell_types <- dat$celltype_090523
clean_cell_types <- str_replace_all(cell_types, "[ /\\-]", "_")
clean_cell_types <- as.factor(clean_cell_types)
names(clean_cell_types) <- colnames(counts)
table(clean_cell_types)
# create emsembl to hgnc dictionary

# marker selection
block <- dat$DataID
marker.info <- scoreMarkers(counts, clean_cell_types, block = block)
saveRDS(marker.info, paste0(save_path, "/blocked_marker.RDS"))
marker.info <- readRDS("/data/xu_lab_projectsx/yuanzhou/drug_resist/results/02_rctd/scran_markers/blocked_marker.RDS")

# create list of selected markers for each cell type
ensg_to_hgnc <- readRDS("/data/xu_lab_projectsx/yuanzhou/drug_resist/data/ensg_to_hgnc_dict.rds")
blocked_marker_list <- list()

# loop through all cell types in marker.info
for (celltype in names(marker.info)) {
  tmp <- marker.info[[celltype]]$mean.AUC
  tmp_genes <- tmp[tmp > 0.75]
  # extract gene names from tmp_genes (which is a named vector)
  blocked_marker_genes <- names(tmp_genes)
  blocked_marker_list[[celltype]] <- blocked_marker_genes
}

# save the list of markers
saveRDS(blocked_marker_list, paste0(save_path, "/blocked_marker_list.RDS"))

# create combined vector of all markers
combined_markers <- unique(unlist(blocked_marker_list))

# remove MT and RB genes
mt_ensg <- readRDS("/data/xu_lab_projectsx/yuanzhou/drug_resist/data/mt_ensg_dict.rds")
rb_ensg <- readRDS("/data/xu_lab_projectsx/yuanzhou/drug_resist/data/rb_ensg_dict.rds")
combined_markers <- combined_markers[!combined_markers %in% mt_ensg]
combined_markers <- combined_markers[!combined_markers %in% rb_ensg]
saveRDS(combined_markers, paste0(save_path, "/blocked_marker_vec.RDS"))

# collect top 2 markers from each cell type for dotplot
blocked_top2_markers <- c()
for (celltype in names(blocked_marker_list)) {
  if (length(blocked_marker_list[[celltype]]) >= 2) {
    # take first 2 markers
    top2_markers <- blocked_marker_list[[celltype]][1:2]
    blocked_top2_markers <- c(blocked_top2_markers, top2_markers)
  } else if (length(blocked_marker_list[[celltype]]) == 1) {
    # if only one marker, use it
    top_marker <- blocked_marker_list[[celltype]][1]
    blocked_top2_markers <- c(blocked_top2_markers, top_marker)
  }
}
# remove duplicates
blocked_top2_markers <- unique(blocked_top2_markers)

# temprory change dat rownames to hgnc symbols
hgnc_to_ensg <- readRDS("/data/xu_lab_projectsx/yuanzhou/drug_resist/data/hgnc_to_ensg_dict.rds")
labels_dict <- hgnc_to_ensg[blocked_top2_markers]
# Remove NA values
labels_dict <- labels_dict[!is.na(labels_dict)]

# Create named vector for scale_x_discrete: ENSG_ID = "HGNC_symbol"
# names are ENSG IDs (what's in the plot), values are HGNC symbols (what to display)
scale_labels <- setNames(names(labels_dict), labels_dict)

# Get ENSG IDs as features
blocked_features_ensg <- as.character(labels_dict)

# do dotplot
dotplot <- DotPlot(dat, 
                   assay = "SoupX", 
                   features = blocked_features_ensg,
                   group.by = "celltype_090523") +
  scale_x_discrete(labels = scale_labels) + 
  theme(panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(save_path, "/blocked_dotplot_all_top2_markers.png"), dotplot, width = 20, height = 8)