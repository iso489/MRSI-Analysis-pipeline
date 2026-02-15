#!/usr/bin/env Rscript
# =============================================================================
# fgsea Gene Set Enrichment Analysis - Fixed for Ensembl ID conversion
# Pre-ranked GSEA against pre-specified metabolic pathways
# =============================================================================

suppressPackageStartupMessages({
  library(fgsea)
  library(data.table)
  library(ggplot2)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

cat("R version:", R.version.string, "\n")
cat("fgsea version:", as.character(packageVersion("fgsea")), "\n\n")

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
BASE_DIR    <- "/scratch/ilyaso/tcga_lgg_open"
DATA_DIR    <- file.path(BASE_DIR, "data_rnaseq_open")
RESULTS_DIR <- file.path(BASE_DIR, "results")

RANKS_FILE <- file.path(DATA_DIR, "ranks.tsv")
GMT_FILE   <- file.path(DATA_DIR, "pathways.gmt")
OUTPUT_FILE <- file.path(RESULTS_DIR, "fgsea_results.tsv")

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

cat("Input ranks:", RANKS_FILE, "\n")
cat("Input GMT  :", GMT_FILE, "\n")
cat("Output file:", OUTPUT_FILE, "\n\n")

# -----------------------------------------------------------------------------
# Step 1: Load ranked gene list
# -----------------------------------------------------------------------------
cat("=== Loading ranked genes ===\n")

if (!file.exists(RANKS_FILE)) {
  stop("Ranks file not found: ", RANKS_FILE,
       "\nRun 01_deseq2_analysis.R first to generate ranks.tsv")
}

ranks_dt <- fread(RANKS_FILE, header = FALSE, col.names = c("gene_id", "stat"))
cat("Loaded", nrow(ranks_dt), "ranked genes\n")

# -----------------------------------------------------------------------------
# Step 2: Convert Ensembl IDs to Gene Symbols
# -----------------------------------------------------------------------------
cat("\n=== Converting Ensembl IDs to Gene Symbols ===\n")

# Remove version numbers if present (ENSG00000188162.10 -> ENSG00000188162)
ranks_dt[, ensembl_id := sub("\\.\\d+$", "", gene_id)]

# Map Ensembl IDs to gene symbols using org.Hs.eg.db
ensembl_ids <- unique(ranks_dt$ensembl_id)
cat("Unique Ensembl IDs:", length(ensembl_ids), "\n")

# Get gene symbols
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = ensembl_ids,
  columns = "SYMBOL",
  keytype = "ENSEMBL"
)
gene_map <- as.data.table(gene_map)
setnames(gene_map, c("ENSEMBL", "SYMBOL"), c("ensembl_id", "symbol"))

# Remove NA mappings and duplicates (keep first)
gene_map <- gene_map[!is.na(symbol)]
gene_map <- gene_map[!duplicated(ensembl_id)]

cat("Mapped genes:", nrow(gene_map), "\n")

# Merge with ranks
ranks_dt <- merge(ranks_dt, gene_map, by = "ensembl_id", all.x = TRUE)
ranks_mapped <- ranks_dt[!is.na(symbol)]

cat("Genes with symbol mapping:", nrow(ranks_mapped), "\n")

# Handle duplicate symbols (keep the one with highest absolute stat)
ranks_mapped[, abs_stat := abs(stat)]
ranks_mapped <- ranks_mapped[order(-abs_stat)]
ranks_mapped <- ranks_mapped[!duplicated(symbol)]
ranks_mapped[, abs_stat := NULL]

cat("Unique gene symbols:", nrow(ranks_mapped), "\n")

# Create named vector for fgsea
ranks <- setNames(ranks_mapped$stat, ranks_mapped$symbol)
ranks <- sort(ranks, decreasing = TRUE)

cat("Rank range:", round(min(ranks), 3), "to", round(max(ranks), 3), "\n")

# Show top genes by symbol
cat("\nTop 5 genes (upregulated in IDH-mutant):\n")
print(head(sort(ranks, decreasing = TRUE), 5))
cat("\nBottom 5 genes (downregulated in IDH-mutant):\n")
print(head(sort(ranks, decreasing = FALSE), 5))

# -----------------------------------------------------------------------------
# Step 3: Load pathway gene sets
# -----------------------------------------------------------------------------
cat("\n=== Loading pathways ===\n")

if (!file.exists(GMT_FILE)) {
  stop("GMT file not found: ", GMT_FILE,
       "\nRun 00_prepare_pathways.R first to create pathways.gmt")
}

pathways <- gmtPathways(GMT_FILE)
cat("Loaded", length(pathways), "pathways:\n")

for (nm in names(pathways)) {
  n_genes <- length(pathways[[nm]])
  n_in_data <- sum(pathways[[nm]] %in% names(ranks))
  cat(sprintf("  %s: %d genes (%d in data)\n", nm, n_genes, n_in_data))
}

# Check if we have gene overlap
total_overlap <- sum(sapply(pathways, function(p) sum(p %in% names(ranks))))
if (total_overlap == 0) {
  stop("No overlap between pathway genes and ranked genes! Check gene ID formats.")
}

# -----------------------------------------------------------------------------
# Step 4: Run fgsea (10,000 permutations as per manuscript)
# -----------------------------------------------------------------------------
cat("\n=== Running fgsea ===\n")

set.seed(42)  # For reproducibility

fgsea_res <- fgsea(
  pathways = pathways,
  stats = ranks,
  minSize = 5,      # Lowered from 10 since some pathways are small
  maxSize = 500,
  nPermSimple = 10000  # As specified in manuscript methods
)

# Sort by p-value
fgsea_res <- fgsea_res[order(pval)]

cat("\nfgsea results:\n")
print(fgsea_res[, .(pathway, pval, padj, NES, size)])

# -----------------------------------------------------------------------------
# Step 5: Interpret results relative to manuscript hypotheses
# -----------------------------------------------------------------------------
cat("\n=== Results interpretation ===\n")
cat("(Positive NES = upregulated in IDH-mutant; Negative NES = downregulated)\n\n")

sig_threshold <- 0.05

for (i in 1:nrow(fgsea_res)) {
  pathway <- fgsea_res$pathway[i]
  nes <- fgsea_res$NES[i]
  padj <- fgsea_res$padj[i]
  size <- fgsea_res$size[i]
  
  direction <- ifelse(nes > 0, "UPREGULATED in IDH-mutant", "DOWNREGULATED in IDH-mutant")
  sig <- ifelse(padj < sig_threshold, "SIGNIFICANT", "not significant")
  
  cat(sprintf("%s:\n", pathway))
  cat(sprintf("  NES = %.3f (%s)\n", nes, direction))
  cat(sprintf("  FDR q = %.4f (%s at q < %.2f)\n", padj, sig, sig_threshold))
  cat(sprintf("  Genes in pathway: %d\n\n", size))
}

# -----------------------------------------------------------------------------
# Step 6: Save results
# -----------------------------------------------------------------------------
cat("\n=== Saving results ===\n")

# Collapse leading edge genes to comma-separated string
fgsea_out <- copy(fgsea_res)
fgsea_out[, leadingEdge := sapply(leadingEdge, paste, collapse = ",")]

fwrite(fgsea_out, OUTPUT_FILE, sep = "\t")
cat("Saved:", OUTPUT_FILE, "\n")

# -----------------------------------------------------------------------------
# Step 7: Generate enrichment plots
# -----------------------------------------------------------------------------
cat("\n=== Generating plots ===\n")

# Enrichment plots for each pathway
for (pathway_name in names(pathways)) {
  if (pathway_name %in% fgsea_res$pathway && fgsea_res[pathway == pathway_name, size] >= 5) {
    pdf_file <- file.path(RESULTS_DIR, paste0("enrichment_", pathway_name, ".pdf"))
    
    p <- plotEnrichment(pathways[[pathway_name]], ranks) +
      labs(title = paste0(pathway_name, "\nNES=", 
                          round(fgsea_res[pathway == pathway_name, NES], 2),
                          ", q=", 
                          format(fgsea_res[pathway == pathway_name, padj], digits = 3))) +
      theme_minimal()
    
    ggsave(pdf_file, p, width = 8, height = 5)
    cat("Saved:", basename(pdf_file), "\n")
  }
}

# Summary bar plot of NES values
if (nrow(fgsea_res) > 0) {
  p_summary <- ggplot(fgsea_res, aes(x = reorder(pathway, NES), y = NES, fill = padj < 0.05)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("grey60", "firebrick"), 
                      labels = c("FDR >= 0.05", "FDR < 0.05"),
                      name = "Significance") +
    labs(title = "Pathway Enrichment: IDH-mutant vs IDH-wildtype",
         subtitle = "Positive NES = upregulated in IDH-mutant",
         x = "Pathway",
         y = "Normalized Enrichment Score (NES)") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  
  ggsave(file.path(RESULTS_DIR, "fgsea_summary.pdf"), p_summary, width = 10, height = 6)
  cat("Saved: fgsea_summary.pdf\n")
}

cat("\n=== fgsea analysis complete ===\n")
