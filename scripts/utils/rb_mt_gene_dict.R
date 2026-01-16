gene_dict <- readRDS("/data/xu_lab_projectsx/yuanzhou/drug_resist/data/ensg_to_hgnc_dict.rds")
hgnc_to_ensg <- readRDS("/data/xu_lab_projectsx/yuanzhou/drug_resist/data/hgnc_to_ensg_dict.rds")
# MT genes
mt_genes <- gene_dict[grep("^MT-", gene_dict)]
rb_genes <- gene_dict[grep("^RP[SL]", gene_dict)]

# 
mt_ensg <- hgnc_to_ensg[mt_genes]
rb_ensg <- hgnc_to_ensg[rb_genes]

# save
saveRDS(mt_ensg, "/data/xu_lab_projectsx/yuanzhou/drug_resist/data/mt_ensg_dict.rds")
saveRDS(rb_ensg, "/data/xu_lab_projectsx/yuanzhou/drug_resist/data/rb_ensg_dict.rds")