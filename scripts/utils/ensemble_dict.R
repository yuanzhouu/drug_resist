library(Seurat)
library(biomaRt)

dat <- readRDS("/data/xu_lab_projectsx/Hasan/VisiumHD/rds/references/lam_lung.rds")

# Get Ensembl IDs from rownames
ensembl_genes <- rownames(dat)

# Connect to Ensembl biomart
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Map Ensembl IDs to HGNC symbols
gene_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), 
                  filters = "ensembl_gene_id", 
                  values = ensembl_genes, 
                  mart = mart)

# Remove rows with empty HGNC symbols
gene_map <- gene_map[gene_map$hgnc_symbol != "", ]

# Create dictionary vector: HGNC -> Ensembl
hgnc_to_ensg <- setNames(gene_map$ensembl_gene_id, gene_map$hgnc_symbol)

# Remove NA mappings
hgnc_to_ensg <- hgnc_to_ensg[!is.na(hgnc_to_ensg)]

# ensg to hgnc dictionary
ensg_to_hgnc <- setNames(gene_map$hgnc_symbol, gene_map$ensembl_gene_id)
ensg_to_hgnc <- ensg_to_hgnc[!is.na(ensg_to_hgnc)]
# Display summary
cat("Dictionary created with", length(hgnc_to_ensg), "HGNC -> Ensembl mappings\n")

# save dictionary
saveRDS(hgnc_to_ensg, "/data/xu_lab_projectsx/yuanzhou/drug_resist/data/hgnc_to_ensg_dict.rds")
saveRDS(ensg_to_hgnc, "/data/xu_lab_projectsx/yuanzhou/drug_resist/data/ensg_to_hgnc_dict.rds")