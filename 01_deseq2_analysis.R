#!/usr/bin/env Rscript
# =============================================================================
# DESeq2 Differential Expression Analysis: IDH-mutant vs IDH-wildtype
# Fixed for GDC STAR gene counts format with UUID subdirectories
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(data.table)
  library(ggplot2)
})

cat("R version:", R.version.string, "\n")
cat("DESeq2 version:", as.character(packageVersion("DESeq2")), "\n\n")

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
BASE_DIR    <- "/scratch/ilyaso/tcga_lgg_open"
DATA_DIR    <- file.path(BASE_DIR, "data_rnaseq_open")
RESULTS_DIR <- file.path(BASE_DIR, "results")
SAMPLE_SHEET <- file.path(BASE_DIR, "tcga_lgg_rnaseq_sample_sheet_annotated.tsv")

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Sample sheet:", SAMPLE_SHEET, "\n")
cat("Data directory:", DATA_DIR, "\n\n")

# -----------------------------------------------------------------------------
# Step 1: Load annotated sample sheet
# -----------------------------------------------------------------------------
cat("=== Loading sample sheet ===\n")

sample_info <- fread(SAMPLE_SHEET)
cat("Total samples in sheet:", nrow(sample_info), "\n")

# Filter to samples with IDH annotation
sample_info <- sample_info[IDH_status %in% c("Mutant", "Wildtype")]
cat("Samples with IDH annotation:", nrow(sample_info), "\n")

cat("\nIDH status distribution:\n")
print(table(sample_info$IDH_status))

# -----------------------------------------------------------------------------
# Step 2: Find count files in UUID subdirectories
# -----------------------------------------------------------------------------
cat("\n=== Finding count files ===\n")

# Find all count files recursively
count_files <- list.files(
  DATA_DIR,
  pattern = "augmented_star_gene_counts\\.tsv$",
  full.names = TRUE,
  recursive = TRUE
)

cat("Found", length(count_files), "count files\n")

# Extract file UUIDs from the file names (the first part before .rna_seq)
# File names look like: a7732b3f-6294-47ec-b869-0aac6ee6c815.rna_seq.augmented_star_gene_counts.tsv
file_uuids <- sub("\\.rna_seq\\.augmented_star_gene_counts\\.tsv$", "", basename(count_files))

# The sample sheet has file_name column with the full filename
# Extract UUID from sample sheet file names
sample_info[, file_uuid := sub("\\.rna_seq\\.augmented_star_gene_counts\\.tsv$", "", count_file)]

# Create mapping from file UUID to full path
file_map <- data.table(
  file_uuid = file_uuids,
  count_path = count_files
)

# Merge to get file paths for each sample
sample_info <- merge(sample_info, file_map, by = "file_uuid", all.x = TRUE)

# Filter to samples with matching count files
sample_info <- sample_info[!is.na(count_path)]
cat("Samples with matching count files:", nrow(sample_info), "\n")

if (nrow(sample_info) < 20) {
  cat("\nDEBUG: Sample file UUIDs (first 5):\n")
  print(head(sample_info$file_uuid, 5))
  cat("\nDEBUG: Count file UUIDs (first 5):\n")
  print(head(file_uuids, 5))
  stop("Too few samples matched. Check UUID matching.")
}

cat("\nFinal sample counts:\n")
print(table(sample_info$IDH_status))

# -----------------------------------------------------------------------------
# Step 3: Read count files and build count matrix
# -----------------------------------------------------------------------------
cat("\n=== Reading count files ===\n")

# Read first file to get gene list
first_file <- sample_info$count_path[1]
cat("Reading first file to get gene structure:", basename(first_file), "\n")

# Skip comment line (starts with #) and read
test_counts <- fread(first_file, skip = 1)  # Skip the "# gene-model" comment line
cat("Columns:", paste(names(test_counts), collapse = ", "), "\n")
cat("Rows:", nrow(test_counts), "\n")

# Use "unstranded" column for total counts
count_col <- "unstranded"
gene_col <- "gene_id"

# Read all count files
cat("\nReading", nrow(sample_info), "count files...\n")

counts_list <- lapply(seq_len(nrow(sample_info)), function(i) {
  f <- sample_info$count_path[i]
  sample_id <- sample_info$patient_id[i]  # Use patient_id as sample identifier
  
  # Read file, skipping comment line
  dt <- fread(f, skip = 1, select = c(gene_col, count_col))
  setnames(dt, c(gene_col, count_col), c("gene_id", sample_id))
  
  if (i %% 25 == 0) cat("  Read", i, "/", nrow(sample_info), "files\n")
  return(dt)
})

cat("  Read", nrow(sample_info), "/", nrow(sample_info), "files\n")

# Merge all counts
cat("\nMerging count matrices...\n")
counts <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), counts_list)

# Convert to matrix
gene_ids <- counts$gene_id
counts_mat <- as.matrix(counts[, -1, with = FALSE])
rownames(counts_mat) <- gene_ids

# Replace NA with 0
counts_mat[is.na(counts_mat)] <- 0

cat("Count matrix dimensions:", nrow(counts_mat), "genes x", ncol(counts_mat), "samples\n")

# -----------------------------------------------------------------------------
# Step 4: Filter genes and prepare for DESeq2
# -----------------------------------------------------------------------------
cat("\n=== Filtering and preparing data ===\n")

# Remove special rows (N_unmapped, N_multimapping, etc.)
special_rows <- grep("^N_", rownames(counts_mat))
if (length(special_rows) > 0) {
  cat("Removing", length(special_rows), "special/summary rows\n")
  counts_mat <- counts_mat[-special_rows, ]
}

# Remove Ensembl version numbers (ENSG00000000003.15 -> ENSG00000000003)
rownames(counts_mat) <- sub("\\.\\d+$", "", rownames(counts_mat))

# Filter low-count genes (keep genes with >= 10 counts in >= 10 samples)
keep <- rowSums(counts_mat >= 10) >= 10
counts_mat <- counts_mat[keep, ]
cat("Genes after filtering:", nrow(counts_mat), "\n")

# Prepare sample metadata (must match column names of counts_mat)
col_data <- data.frame(
  sample_id = colnames(counts_mat),
  row.names = colnames(counts_mat)
)

# Match IDH status to samples
col_data$IDH_status <- sample_info$IDH_status[match(col_data$sample_id, sample_info$patient_id)]

# Add grade if available
if ("WHO_grade" %in% names(sample_info)) {
  col_data$WHO_grade <- sample_info$WHO_grade[match(col_data$sample_id, sample_info$patient_id)]
}

# Remove samples without IDH annotation (shouldn't happen, but just in case)
valid_samples <- !is.na(col_data$IDH_status)
counts_mat <- counts_mat[, valid_samples]
col_data <- col_data[valid_samples, , drop = FALSE]

cat("\nFinal dataset:\n")
cat("  Samples:", nrow(col_data), "\n")
cat("  IDH-Mutant:", sum(col_data$IDH_status == "Mutant"), "\n")
cat("  IDH-Wildtype:", sum(col_data$IDH_status == "Wildtype"), "\n")

# Set factor levels (Wildtype as reference)
col_data$IDH_status <- factor(col_data$IDH_status, levels = c("Wildtype", "Mutant"))

# -----------------------------------------------------------------------------
# Step 5: Run DESeq2
# -----------------------------------------------------------------------------
cat("\n=== Running DESeq2 ===\n")

# Build design formula
design_formula <- ~ IDH_status

# Add grade as covariate if available and has variation
if ("WHO_grade" %in% names(col_data)) {
  grade_valid <- !is.na(col_data$WHO_grade)
  if (sum(grade_valid) > 0 && length(unique(col_data$WHO_grade[grade_valid])) > 1) {
    col_data$WHO_grade <- factor(col_data$WHO_grade)
    design_formula <- ~ WHO_grade + IDH_status
    cat("Design: ~ WHO_grade + IDH_status\n")
  } else {
    cat("Design: ~ IDH_status (grade not usable)\n")
  }
} else {
  cat("Design: ~ IDH_status\n")
}

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_mat),  # Ensure integers
  colData = col_data,
  design = design_formula
)

# Run DESeq2
cat("Running DESeq2 (this may take a few minutes)...\n")
dds <- DESeq(dds)

# Extract results
res <- results(dds, contrast = c("IDH_status", "Mutant", "Wildtype"))
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

cat("\nDESeq2 results summary:\n")
summary(res)

# -----------------------------------------------------------------------------
# Step 6: Generate ranked gene list for fgsea
# -----------------------------------------------------------------------------
cat("\n=== Generating ranked gene list for fgsea ===\n")

# Rank by Wald statistic (as specified in manuscript methods)
ranks_df <- data.table(
  gene_id = res_df$gene_id,
  stat = res_df$stat
)

# Remove NA values
ranks_df <- ranks_df[!is.na(stat)]

# Sort by statistic (descending)
ranks_df <- ranks_df[order(-stat)]

cat("Ranked genes:", nrow(ranks_df), "\n")
cat("\nTop 5 genes (upregulated in IDH-mutant):\n")
print(head(ranks_df, 5))
cat("\nBottom 5 genes (downregulated in IDH-mutant):\n")
print(tail(ranks_df, 5))

# Save ranks file - THIS IS WHAT fgsea NEEDS
ranks_file <- file.path(DATA_DIR, "ranks.tsv")
fwrite(ranks_df, ranks_file, sep = "\t", col.names = FALSE)
cat("\n*** Saved ranks file:", ranks_file, "***\n")

# Verify the file was created
if (file.exists(ranks_file)) {
  cat("Verified: ranks.tsv exists with", nrow(fread(ranks_file)), "genes\n")
} else {
  stop("ERROR: Failed to create ranks.tsv!")
}

# -----------------------------------------------------------------------------
# Step 7: Save additional outputs
# -----------------------------------------------------------------------------
cat("\n=== Saving outputs ===\n")

# Full DESeq2 results
results_file <- file.path(RESULTS_DIR, "deseq2_results.tsv")
fwrite(res_df, results_file, sep = "\t")
cat("Saved:", results_file, "\n")

# VST-normalized counts for GSVA
cat("Computing VST normalization...\n")
vst <- varianceStabilizingTransformation(dds, blind = FALSE)
saveRDS(assay(vst), file.path(RESULTS_DIR, "vst_counts.rds"))
cat("Saved: vst_counts.rds\n")

# Sample metadata
saveRDS(col_data, file.path(RESULTS_DIR, "sample_metadata.rds"))
cat("Saved: sample_metadata.rds\n")

# Volcano plot
cat("Generating volcano plot...\n")
pdf(file.path(RESULTS_DIR, "volcano_plot.pdf"), width = 8, height = 6)
volcano_df <- res_df[!is.na(res_df$padj), ]
volcano_df$significant <- volcano_df$padj < 0.05 & abs(volcano_df$log2FoldChange) > 1
plot(volcano_df$log2FoldChange, -log10(volcano_df$pvalue),
     col = ifelse(volcano_df$significant, "firebrick", "grey60"),
     pch = 20, cex = 0.5,
     xlab = "log2 Fold Change (Mutant vs Wildtype)",
     ylab = "-log10(p-value)",
     main = "Differential Expression: IDH-mutant vs IDH-wildtype")
abline(h = -log10(0.05), lty = 2, col = "grey40")
abline(v = c(-1, 1), lty = 2, col = "grey40")
n_up <- sum(volcano_df$significant & volcano_df$log2FoldChange > 0, na.rm = TRUE)
n_down <- sum(volcano_df$significant & volcano_df$log2FoldChange < 0, na.rm = TRUE)
legend("topright", legend = c(paste("Up in Mutant:", n_up), paste("Down in Mutant:", n_down)),
       pch = 20, col = "firebrick", bty = "n")
dev.off()
cat("Saved: volcano_plot.pdf\n")

cat("\n=== DESeq2 analysis complete ===\n")
cat("Output files in:", RESULTS_DIR, "\n")
cat("Ranks file ready for fgsea:", ranks_file, "\n")
