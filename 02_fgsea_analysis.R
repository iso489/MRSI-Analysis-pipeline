!/usr/bin/env Rscript
# Script 02: fgsea — Pre-ranked GSEA, glioma vs healthy normal brain
#
# INPUTS  (from script 01 and 00)
#   BASE_DIR/results/glioma_vs_normal_ranked_genes.rds  — named vector of
#     DESeq2 Wald statistics, ENSG IDs as names, sorted descending
#     stat > 0 = higher in glioma; stat < 0 = lower in glioma (higher in normal)
#   BASE_DIR/metadata/ensg_to_symbol.rds  — ENSG → HGNC symbol map
#   BASE_DIR/metadata/pathways.gmt        — from script 00

suppressPackageStartupMessages({
  library(fgsea)
  library(data.table)
  library(ggplot2)
})

cat("fgsea version:", as.character(packageVersion("fgsea")), "\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================
USER     <- Sys.getenv("USER")
BASE_DIR <- file.path("/scratch", USER, "glioma_vs_normal")
RES_DIR  <- file.path(BASE_DIR, "results")
META_DIR <- file.path(BASE_DIR, "metadata")

RANKED_FILE <- file.path(RES_DIR,  "glioma_vs_normal_ranked_genes.rds")
ENSG_MAP    <- file.path(META_DIR, "ensg_to_symbol.rds")
GMT_FILE    <- file.path(META_DIR, "pathways.gmt")
OUTPUT_FILE <- file.path(RES_DIR,  "fgsea_results.tsv")

for (f in c(RANKED_FILE, ENSG_MAP, GMT_FILE))
  if (!file.exists(f)) stop(sprintf("Required file not found: %s", f))

# =============================================================================
# STEP 1: Load ranked gene list and convert ENSG → symbol
#
# Script 01 saves a named numeric vector (names = ENSG IDs, values = Wald stat).
# We convert to symbols here using the SAME map script 01 built, guaranteeing
# the symbols are consistent with the pathway gene sets in the GMT file.
# =============================================================================
cat("--- Step 1: Loading ranked genes ---\n")

ranked_ensg <- readRDS(RANKED_FILE)
cat(sprintf("  Loaded: %d ranked genes\n", length(ranked_ensg)))
cat(sprintf("  Stat range: [%.3f, %.3f]\n", min(ranked_ensg), max(ranked_ensg)))

ensg_map <- readRDS(ENSG_MAP)
cat(sprintf("  ENSG map: %d entries\n", nrow(ensg_map)))

# Map ENSG → symbol
id_vec    <- setNames(ensg_map$symbol, ensg_map$gene_id)
symbols   <- id_vec[names(ranked_ensg)]
mapped    <- !is.na(symbols)
cat(sprintf("  Genes with symbol: %d / %d\n", sum(mapped), length(ranked_ensg)))

ranked_sym <- ranked_ensg[mapped]
names(ranked_sym) <- symbols[mapped]

# If multiple ENSG IDs map to same symbol, keep the one with highest |stat|
dt <- data.table(symbol = names(ranked_sym), stat = ranked_sym)
dt <- dt[order(-abs(stat))]
dt <- dt[!duplicated(symbol)]
ranked <- setNames(dt$stat, dt$symbol)
ranked <- sort(ranked, decreasing = TRUE)

cat(sprintf("  Unique symbols after dedup: %d\n\n", length(ranked)))
cat("  Top 5 (upregulated in glioma):\n")
print(head(ranked, 5))
cat("  Bottom 5 (downregulated in glioma / higher in normal brain):\n")
print(tail(ranked, 5))

# =============================================================================
# STEP 2: Load pathways
# =============================================================================
cat("\n--- Step 2: Loading pathways ---\n")

pathways <- gmtPathways(GMT_FILE)
cat(sprintf("  Loaded %d pathways:\n", length(pathways)))

for (nm in names(pathways)) {
  n_total <- length(pathways[[nm]])
  n_in    <- sum(pathways[[nm]] %in% names(ranked))
  cat(sprintf("    %-50s %d genes (%d in data)\n", nm, n_total, n_in))
}

if (sum(sapply(pathways, function(p) sum(p %in% names(ranked)))) == 0L)
  stop("No overlap between pathway genes and ranked list — check symbol mapping.")

# =============================================================================
# STEP 3: Run fgsea
# =============================================================================
cat("\n--- Step 3: Running fgsea (10,000 permutations) ---\n")

set.seed(42)
fgsea_res <- fgsea(
  pathways    = pathways,
  stats       = ranked,
  minSize     = 5,
  maxSize     = 500,
  nPermSimple = 10000
)
fgsea_res <- fgsea_res[order(pval)]

cat("\n  Results:\n")
print(fgsea_res[, .(pathway, pval, padj, NES, size)])

# =============================================================================
# STEP 4: Interpret results
# =============================================================================
cat("\n--- Step 4: Interpretation ---\n")
cat("  NES > 0 = upregulated in glioma vs normal brain\n")
cat("  NES < 0 = downregulated in glioma (higher in normal brain)\n\n")

for (i in seq_len(nrow(fgsea_res))) {
  pw   <- fgsea_res$pathway[i]
  nes  <- fgsea_res$NES[i]
  q    <- fgsea_res$padj[i]
  sz   <- fgsea_res$size[i]
  dir  <- ifelse(nes > 0, "UPREGULATED in glioma", "DOWNREGULATED in glioma")
  sig  <- ifelse(q < 0.05, "SIGNIFICANT", "not significant")
  cat(sprintf("%s:\n  NES = %.3f (%s)\n  FDR q = %.4f (%s)\n  Genes: %d\n\n",
              pw, nes, dir, q, sig, sz))
}

# =============================================================================
# STEP 5: Save results
# =============================================================================
cat("--- Step 5: Saving results ---\n")

fgsea_out <- copy(fgsea_res)
fgsea_out[, leadingEdge := sapply(leadingEdge, paste, collapse = ",")]
fwrite(fgsea_out, OUTPUT_FILE, sep = "\t")
cat(sprintf("  Saved: %s\n", basename(OUTPUT_FILE)))

# =============================================================================
# STEP 6: Plots
# =============================================================================
cat("\n--- Step 6: Generating plots ---\n")

# Per-pathway enrichment plots
for (pw in names(pathways)) {
  row <- fgsea_res[pathway == pw]
  if (nrow(row) == 0L || row$size < 5L) next
  pdf_file <- file.path(RES_DIR, sprintf("enrichment_%s.pdf", pw))
  p <- plotEnrichment(pathways[[pw]], ranked) +
    labs(title = pw,
         subtitle = sprintf("NES = %.2f  |  FDR q = %.4f  |  n = %d genes",
                            row$NES, row$padj, row$size)) +
    theme_minimal()
  ggsave(pdf_file, p, width = 8, height = 5)
  cat(sprintf("  Saved: %s\n", basename(pdf_file)))
}

# Summary bar plot
if (nrow(fgsea_res) > 0L) {
  p_sum <- ggplot(fgsea_res,
                  aes(x = reorder(pathway, NES), y = NES,
                      fill = padj < 0.05)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("grey60", "firebrick"),
                      labels = c("FDR \u2265 0.05", "FDR < 0.05"),
                      name   = "Significance") +
    labs(title    = "Pathway Enrichment: Glioma vs Normal Brain",
         subtitle = "NES > 0 = upregulated in glioma",
         x = NULL, y = "Normalized Enrichment Score (NES)") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))
  ggsave(file.path(RES_DIR, "fgsea_summary.pdf"), p_sum,
         width = 10, height = 6)
  cat("  Saved: fgsea_summary.pdf\n")
}

cat("\n=== Script 02 complete ===\n")
cat(sprintf("Results: %s\n", RES_DIR))# -----------------------------------------------------------------------------
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
