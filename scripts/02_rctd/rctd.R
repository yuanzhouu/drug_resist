library(Seurat)
library(spacexr)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide sample name as argument")
}
sample <- args[1]

# Define paths
ref_path <- "/data/xu_lab_projectsx/yuanzhou/drug_resist/results/02_rctd/ref/ref_lam_lung.rds"
query_path <- paste0("/data/xu_lab_projectsx/yuanzhou/drug_resist/results/02_rctd/query/query_", sample, ".rds")
output_dir <- "/data/xu_lab_projectsx/yuanzhou/drug_resist/results/02_rctd/test_results"
sample_dir <- paste0(output_dir, "/", sample)
dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)

# Redirect all output to log file
log_file <- paste0(sample_dir, "/RCTD_", sample, "_log.txt")
sink(log_file, split = TRUE)  # split = TRUE keeps output in console and file

cat("========================================\n")
cat("RCTD Analysis for sample:", sample, "\n")
cat("Start time:", as.character(Sys.time()), "\n")
cat("========================================\n\n")

# Load reference and query objects
cat("Loading reference object...\n")
cat("Reference path:", ref_path, "\n")
reference <- readRDS(ref_path)
cat("Reference loaded successfully.\n\n")

cat("Loading query object for sample:", sample, "\n")
cat("Query path:", query_path, "\n")
query <- readRDS(query_path)
cat("Query loaded successfully.\n\n")

# Run RCTD
cat("Creating RCTD object...\n")
cat("Parameters: UMI_min = 10, max_cores = 20\n")
myRCTD <- create.RCTD(query, reference, UMI_min = 10, UMI_min_sigma = 100, max_cores = 20)
cat("RCTD object created successfully.\n\n")

cat("Running RCTD with doublet mode...\n")
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'doublet')
cat("RCTD analysis completed.\n\n")

# Save results
output_file <- paste0(sample_dir, "/RCTD_doublet_", sample, ".tsv")
cat("Saving results table to:", output_file, "\n")
write.table(
    myRCTD@results$results_df,
    output_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)
cat("Results table saved.\n\n")

# Save full RCTD object
rctd_output_file <- paste0(sample_dir, "/RCTD_", sample, ".rds")
cat("Saving full RCTD object to:", rctd_output_file, "\n")
saveRDS(myRCTD, file = rctd_output_file)
cat("RCTD object saved.\n\n")

cat("========================================\n")
cat("Completed processing for sample:", sample, "\n")
cat("End time:", as.character(Sys.time()), "\n")
cat("========================================\n")

# Stop redirecting output
sink()