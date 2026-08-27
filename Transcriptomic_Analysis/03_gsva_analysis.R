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
