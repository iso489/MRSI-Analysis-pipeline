#!/usr/bin/env Rscript
# Script 03: GSVA — Per-sample pathway activity, glioma vs normal brain
#
# INPUTS  (from scripts 00 and 01)
#   BASE_DIR/results/DESeq2_dds_glioma_vs_normal.rds  — DESeqDataSet
#   BASE_DIR/metadata/ensg_to_symbol.rds              — ENSG → symbol map
#   BASE_DIR/metadata/pathways.gmt                    — from script 00
#
# METHOD
#   VST normalization is computed from the saved DESeqDataSet (blind = FALSE,
#   preserving the condition structure). ENSG IDs → symbols via the same map
#   used in scripts 01 and 02, guaranteeing consistency. Duplicate symbols are
#   averaged (mean VST) rather than dropped.

suppressPackageStartupMessages({
  library(DESeq2)
  library(GSVA)
  library(data.table)
})

# Auto-install pheatmap / RColorBrewer if needed
for (pkg in c("pheatmap", "RColorBrewer")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cran.mirror.rafal.ca")
}
suppressPackageStartupMessages({
  library(pheatmap)
  library(RColorBrewer)
})

cat("GSVA version:", as.character(packageVersion("GSVA")), "\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================
USER     <- Sys.getenv("USER")
BASE_DIR <- file.path("/scratch", USER, "glioma_vs_normal")
RES_DIR  <- file.path(BASE_DIR, "results")
META_DIR <- file.path(BASE_DIR, "metadata")

DDS_FILE  <- file.path(RES_DIR,  "DESeq2_dds_glioma_vs_normal.rds")
ENSG_MAP  <- file.path(META_DIR, "ensg_to_symbol.rds")
GMT_FILE  <- file.path(META_DIR, "pathways.gmt")

for (f in c(DDS_FILE, ENSG_MAP, GMT_FILE))
  if (!file.exists(f)) stop(sprintf("Required file not found: %s", f))

# =============================================================================
# STEP 1: VST normalization from saved DESeqDataSet
# =============================================================================
cat("--- Step 1: VST normalization ---\n")

dds  <- readRDS(DDS_FILE)
cat(sprintf("  Loaded DDS: %d genes × %d samples\n",
            nrow(dds), ncol(dds)))

vst      <- varianceStabilizingTransformation(dds, blind = FALSE)
expr_mat <- assay(vst)
cat(sprintf("  VST matrix: %d genes × %d samples\n",
            nrow(expr_mat), ncol(expr_mat)))

# Extract condition from colData (glioma / normal)
condition <- as.character(dds$condition)
names(condition) <- colnames(dds)

# =============================================================================
# STEP 2: ENSG → symbol conversion
#
# Uses the same map as scripts 01 and 02 — symbols are guaranteed consistent
# with the pathway GMT file.
# =============================================================================
cat("\n--- Step 2: ENSG → symbol conversion ---\n")

ensg_map <- readRDS(ENSG_MAP)
id_vec   <- setNames(ensg_map$symbol, ensg_map$gene_id)

# Strip version suffix if present (shouldn't be after script 01, but defensive)
row_ids  <- sub("\\.[0-9]+$", "", rownames(expr_mat))
symbols  <- id_vec[row_ids]
keep     <- !is.na(symbols)
cat(sprintf("  Genes with symbol: %d / %d\n", sum(keep), nrow(expr_mat)))

expr_keep    <- expr_mat[keep, , drop = FALSE]
syms_keep    <- symbols[keep]

# Average rows that map to the same symbol
dt     <- as.data.table(expr_keep)
dt[, symbol := syms_keep]
dt_agg <- dt[, lapply(.SD, mean), by = symbol, .SDcols = colnames(expr_keep)]
expr_sym              <- as.matrix(dt_agg[, -"symbol"])
rownames(expr_sym)    <- dt_agg$symbol
cat(sprintf("  Unique symbols after averaging: %d\n", nrow(expr_sym)))

# =============================================================================
# STEP 3: Load pathways
# =============================================================================
cat("\n--- Step 3: Loading pathways ---\n")

gmt_lines <- readLines(GMT_FILE)
pathways  <- lapply(gmt_lines, function(l) {
  p <- strsplit(l, "\t")[[1]]
  p[3:length(p)]
})
names(pathways) <- vapply(gmt_lines, function(l) strsplit(l, "\t")[[1]][1],
                          character(1))

pathways_filt <- lapply(pathways,
                        function(g) intersect(g, rownames(expr_sym)))
cat(sprintf("  Loaded %d pathways:\n", length(pathways_filt)))
for (nm in names(pathways_filt))
  cat(sprintf("    %-50s %d genes in data\n", nm,
              length(pathways_filt[[nm]])))

pathways_filt <- pathways_filt[sapply(pathways_filt, length) >= 5L]
cat(sprintf("  Pathways with >= 5 genes: %d\n", length(pathways_filt)))
if (length(pathways_filt) == 0L)
  stop("No pathways have sufficient genes in the expression data.")

# =============================================================================
# STEP 4: Run GSVA
# =============================================================================
cat("\n--- Step 4: Running GSVA (Gaussian kernel) ---\n")

gsva_param  <- gsvaParam(exprData  = expr_sym,
                          geneSets  = pathways_filt,
                          kcdf      = "Gaussian")
gsva_scores <- gsva(gsva_param, verbose = TRUE)
cat(sprintf("\n  GSVA scores: %d pathways × %d samples\n",
            nrow(gsva_scores), ncol(gsva_scores)))

# =============================================================================
# STEP 5: Statistical comparison — glioma vs normal (Wilcoxon)
#
# Using Wilcoxon rank-sum rather than t-test: GSVA scores are bounded [-1, 1]
# and not guaranteed normal. Wilcoxon is more appropriate for this distribution.
# =============================================================================
cat("\n--- Step 5: Glioma vs normal comparison (Wilcoxon) ---\n")

glioma_ids <- names(condition)[condition == "glioma"]
normal_ids <- names(condition)[condition == "normal"]
cat(sprintf("  Glioma: %d | Normal: %d\n",
            length(glioma_ids), length(normal_ids)))

results_list <- lapply(rownames(gsva_scores), function(pw) {
  sc       <- gsva_scores[pw, ]
  g_scores <- sc[glioma_ids]
  n_scores <- sc[normal_ids]
  wt       <- wilcox.test(g_scores, n_scores, exact = FALSE)
  data.table(
    pathway        = pw,
    mean_glioma    = mean(g_scores),
    mean_normal    = mean(n_scores),
    diff           = mean(g_scores) - mean(n_scores),
    W_stat         = as.numeric(wt$statistic),
    pvalue         = as.numeric(wt$p.value)
  )
})

gsva_results       <- rbindlist(results_list)
gsva_results[, padj := p.adjust(pvalue, method = "BH")]
gsva_results        <- gsva_results[order(pvalue)]

cat("\n  Pathway activity — glioma vs normal:\n")
print(gsva_results[, .(pathway, diff, pvalue, padj)])

fwrite(gsva_results, file.path(RES_DIR, "gsva_statistics.tsv"), sep = "\t")

gsva_out <- as.data.table(t(gsva_scores), keep.rownames = "sample_id")
fwrite(gsva_out, file.path(RES_DIR, "gsva_scores.tsv"), sep = "\t")
cat("\n  Saved: gsva_statistics.tsv | gsva_scores.tsv\n")

# =============================================================================
# STEP 6: Heatmap
# =============================================================================
cat("\n--- Step 6: Heatmap ---\n")

# Order samples: normal first, glioma second
sample_order <- c(names(condition)[condition == "normal"],
                  names(condition)[condition == "glioma"])
annotation_col <- data.frame(
  Condition = condition[sample_order],
  row.names = sample_order
)
ann_colors <- list(Condition = c(normal = "#377EB8", glioma = "#E41A1C"))

pdf(file.path(RES_DIR, "gsva_heatmap.pdf"), width = 14, height = 6)
pheatmap(
  gsva_scores[, sample_order, drop = FALSE],
  annotation_col  = annotation_col,
  annotation_colors = ann_colors,
  cluster_cols    = FALSE,
  cluster_rows    = TRUE,
  scale           = "none",
  fontsize_row    = 10,
  fontsize_col    = 3,
  show_colnames   = FALSE,
  main            = "GSVA Pathway Activity: Glioma vs Normal Brain"
)
dev.off()
cat("  Saved: gsva_heatmap.pdf\n")

# =============================================================================
# STEP 7: Boxplots
# =============================================================================
cat("\n--- Step 7: Boxplots ---\n")

n_pw   <- nrow(gsva_scores)
n_cols <- min(3L, n_pw)
n_rows <- ceiling(n_pw / n_cols)

pdf(file.path(RES_DIR, "gsva_boxplots.pdf"),
    width = 5 * n_cols, height = 5 * n_rows)
par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 4, 1))

for (pw in rownames(gsva_scores)) {
  scores_df <- data.frame(
    score     = as.numeric(gsva_scores[pw, sample_order]),
    Condition = condition[sample_order]
  )
  q    <- gsva_results[pathway == pw, padj][1L]
  praw <- gsva_results[pathway == pw, pvalue][1L]

  boxplot(score ~ Condition, data = scores_df,
          main = gsub("_", "\n", pw),
          ylab = "GSVA Score",
          col  = c(normal = "#377EB8", glioma = "#E41A1C"),
          cex.main = 0.75)
  mtext(sprintf("p = %.3g  |  FDR q = %.3g", praw, q),
        side = 3, line = 0.2, cex = 0.65)
}
dev.off()
cat("  Saved: gsva_boxplots.pdf\n")

cat("\n=== Script 03 complete ===\n")
cat(sprintf("Results: %s\n", RES_DIR))# -----------------------------------------------------------------------------
BASE_DIR    <- "/scratch/ilyaso/tcga_lgg_open"
DATA_DIR    <- file.path(BASE_DIR, "data_rnaseq_open")
RESULTS_DIR <- file.path(BASE_DIR, "results")
GMT_FILE    <- file.path(DATA_DIR, "pathways.gmt")

dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# Step 1: Load normalized expression data
# -----------------------------------------------------------------------------
cat("=== Loading expression data ===\n")

vst_file <- file.path(RESULTS_DIR, "vst_counts.rds")
if (!file.exists(vst_file)) stop("VST counts not found. Run 01_deseq2_analysis.R first.")

expr_mat <- readRDS(vst_file)
cat("Expression matrix:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n")

# -----------------------------------------------------------------------------
# Step 2: Convert Ensembl IDs to Gene Symbols (and average duplicates)
# -----------------------------------------------------------------------------
cat("\n=== Converting Ensembl IDs to Gene Symbols ===\n")

ensembl_ids <- sub("\\.\\d+$", "", rownames(expr_mat))

gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(ensembl_ids),
  columns = "SYMBOL",
  keytype = "ENSEMBL"
)
gene_map <- as.data.table(gene_map)
setnames(gene_map, c("ENSEMBL", "SYMBOL"), c("ensembl_id", "symbol"))
gene_map <- gene_map[!is.na(symbol)]
gene_map <- gene_map[!duplicated(ensembl_id)]

id_to_symbol <- setNames(gene_map$symbol, gene_map$ensembl_id)
symbols <- id_to_symbol[ensembl_ids]
keep <- !is.na(symbols)

expr_keep <- expr_mat[keep, , drop = FALSE]
symbols_keep <- symbols[keep]

# Average duplicate symbols properly
dt <- as.data.table(expr_keep)
dt[, symbol := symbols_keep]
dt_agg <- dt[, lapply(.SD, mean), by = symbol]
expr_mat_symbol <- as.matrix(dt_agg[, -"symbol"])
rownames(expr_mat_symbol) <- dt_agg$symbol

cat("Genes with symbol mapping (after averaging duplicates):", nrow(expr_mat_symbol), "\n")

# -----------------------------------------------------------------------------
# Step 3: Load sample metadata
# -----------------------------------------------------------------------------
cat("\n=== Loading sample metadata ===\n")

metadata_file <- file.path(RESULTS_DIR, "sample_metadata.rds")
if (file.exists(metadata_file)) {
  sample_info <- readRDS(metadata_file)
  sample_info <- as.data.table(sample_info, keep.rownames = "sample_id")
} else {
  sample_sheet <- file.path(BASE_DIR, "tcga_lgg_rnaseq_sample_sheet_annotated.tsv")
  sample_info <- fread(sample_sheet)
  sample_info <- sample_info[IDH_status %in% c("Mutant", "Wildtype")]
}

cat("Samples in metadata:", nrow(sample_info), "\n")

common_samples <- intersect(colnames(expr_mat_symbol), sample_info$sample_id)
cat("Matched samples:", length(common_samples), "\n")
if (length(common_samples) == 0) stop("No samples match between expression matrix and metadata.")

expr_mat_symbol <- expr_mat_symbol[, common_samples, drop = FALSE]
sample_info <- sample_info[sample_id %in% common_samples]
sample_info <- sample_info[match(common_samples, sample_id)]

cat("\nIDH status distribution:\n")
print(table(sample_info$IDH_status))

# -----------------------------------------------------------------------------
# Step 4: Load pathway gene sets
# -----------------------------------------------------------------------------
cat("\n=== Loading pathways ===\n")

if (!file.exists(GMT_FILE)) stop("GMT file not found. Run 00_prepare_pathways.R first.")

gmt_lines <- readLines(GMT_FILE)
pathways <- lapply(gmt_lines, function(line) {
  parts <- strsplit(line, "\t")[[1]]
  genes <- parts[3:length(parts)]
  genes[genes != ""]
})
names(pathways) <- vapply(gmt_lines, function(line) strsplit(line, "\t")[[1]][1], character(1))

cat("Loaded", length(pathways), "pathways\n")

pathways_filtered <- lapply(pathways, function(genes) intersect(genes, rownames(expr_mat_symbol)))

for (nm in names(pathways_filtered)) {
  cat(sprintf("  %s: %d genes in data\n", nm, length(pathways_filtered[[nm]])))
}

pathways_filtered <- pathways_filtered[sapply(pathways_filtered, length) >= 5]
cat("\nPathways with >= 5 genes:", length(pathways_filtered), "\n")
if (length(pathways_filtered) == 0) stop("No pathways have sufficient genes in the expression data.")

# -----------------------------------------------------------------------------
# Step 5: Run GSVA
# -----------------------------------------------------------------------------
cat("\n=== Running GSVA ===\n")

gsva_param <- gsvaParam(
  exprData = expr_mat_symbol,
  geneSets = pathways_filtered,
  kcdf = "Gaussian"
)

gsva_scores <- gsva(gsva_param, verbose = TRUE)
cat("\nGSVA scores computed:", nrow(gsva_scores), "pathways x", ncol(gsva_scores), "samples\n")

# -----------------------------------------------------------------------------
# Step 6: Statistical comparison between IDH groups
# -----------------------------------------------------------------------------
cat("\n=== Comparing pathway activity by IDH status ===\n")

mutant_ids <- sample_info[IDH_status == "Mutant", sample_id]
wild_ids   <- sample_info[IDH_status == "Wildtype", sample_id]

results_list <- vector("list", length = nrow(gsva_scores))
names(results_list) <- rownames(gsva_scores)

for (pathway in rownames(gsva_scores)) {
  scores <- gsva_scores[pathway, ]

  mutant_scores <- scores[mutant_ids]
  wild_scores   <- scores[wild_ids]

  test_res <- t.test(mutant_scores, wild_scores)

  results_list[[pathway]] <- data.table(
    pathway = pathway,
    mean_mutant = mean(mutant_scores),
    mean_wildtype = mean(wild_scores),
    diff = mean(mutant_scores) - mean(wild_scores),
    t_stat = as.numeric(test_res$statistic),
    pvalue = as.numeric(test_res$p.value),
    ci_lower = test_res$conf.int[1],
    ci_upper = test_res$conf.int[2]
  )
}

gsva_results <- rbindlist(results_list, use.names = TRUE, fill = TRUE)
gsva_results[, padj := p.adjust(pvalue, method = "BH")]
gsva_results <- gsva_results[order(pvalue)]

cat("\nGSVA pathway activity differences:\n")
print(gsva_results[, .(pathway, diff, pvalue, padj)])

fwrite(gsva_results, file.path(RESULTS_DIR, "gsva_statistics.tsv"), sep = "\t")
cat("\nSaved: gsva_statistics.tsv\n")

gsva_out <- as.data.table(t(gsva_scores), keep.rownames = "sample_id")
fwrite(gsva_out, file.path(RESULTS_DIR, "gsva_scores.tsv"), sep = "\t")
cat("Saved: gsva_scores.tsv\n")

# -----------------------------------------------------------------------------
# Step 7: Heatmap + boxplots
# -----------------------------------------------------------------------------
cat("\n=== Generating heatmap ===\n")

annotation_col <- data.frame(
  IDH_status = sample_info$IDH_status,
  row.names = sample_info$sample_id
)

# Order samples by IDH status
sample_order <- sample_info[order(IDH_status), sample_id]

pdf(file.path(RESULTS_DIR, "gsva_heatmap.pdf"), width = 12, height = 6)
pheatmap(
  gsva_scores[, sample_order, drop = FALSE],
  annotation_col = annotation_col,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  scale = "none",
  fontsize_row = 10,
  fontsize_col = 4,
  main = "GSVA Pathway Activity Scores by IDH Status",
  show_colnames = FALSE
)
dev.off()
cat("Saved: gsva_heatmap.pdf\n")

cat("\n=== Generating boxplots ===\n")

# Boxplot comparison
pdf(file.path(RESULTS_DIR, "gsva_boxplots.pdf"), width = 12, height = 8)

n_pathways <- nrow(gsva_scores)
n_cols <- min(3, n_pathways)
n_rows <- ceiling(n_pathways / n_cols)
par(mfrow = c(n_rows, n_cols), mar = c(4, 4, 3, 1))

for (pw in rownames(gsva_scores)) {
  scores_df <- data.frame(
    score = as.numeric(gsva_scores[pw, sample_info$sample_id]),
    IDH   = sample_info$IDH_status
  )

  pval <- gsva_results[pathway == pw, pvalue][1]

  boxplot(score ~ IDH, data = scores_df,
          main = gsub("_", "\n", pw),
          ylab = "GSVA Score",
          col = c("#E41A1C", "#377EB8"),
          cex.main = 0.8)

  mtext(sprintf("p = %.3g", pval), side = 3, line = 0, cex = 0.7)
}

dev.off()
cat("Saved: gsva_boxplots.pdf\n")

cat("\n=== GSVA analysis complete ===\n")
cat("Results saved to:", RESULTS_DIR, "\n")
