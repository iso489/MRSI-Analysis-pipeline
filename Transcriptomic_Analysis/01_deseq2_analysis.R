#!/usr/bin/env Rscript
# Script 01: DESeq2 — Glioma vs Healthy Normal Brain
#
# DATA SOURCE
#   UCSC Xena / Toil re-compute hub (Vivian et al., Nat. Biotechnol. 2017)
#   All TCGA + GTEx samples processed with identical STAR/RSEM pipeline.
#
# REQUIRED INPUT FILES (already downloaded):
#   /scratch/ilyaso/glioma_vs_normal/TcgaTargetGtex_gene_expected_count.gz
#   /scratch/ilyaso/glioma_vs_normal/TcgaTargetGTEX_phenotype.txt.gz
#
# NOTE: Matrix stores log2(RSEM_expected_count + 1) -> back-transformed here.
#
# SAMPLE SELECTION
#   TCGA : Primary Tumor only; TCGA-GBM and TCGA-LGG  (condition = "glioma")
#   GTEx : THREE brain regions combined (condition = "normal"):
#            - Brain - Frontal Cortex (BA9)
#            - Brain - Anterior Cingulate Cortex (BA24)
#            - Brain - Cortex
#          Matches the GCT files originally downloaded:
#            gene_reads_v11_brain_frontal_cortex_ba9.gct
#            gene_reads_v11_brain_anterior_cingulate_cortex_ba24.gct
#            gene_reads_v11_brain_cortex.gct

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(BiocParallel)
})

# =============================================================================
# CONFIGURATION
# =============================================================================
BASE_DIR   <- "/scratch/ilyaso/glioma_vs_normal"
RES_DIR    <- file.path(BASE_DIR, "results")
META_DIR   <- file.path(BASE_DIR, "metadata")
COUNT_FILE <- file.path(BASE_DIR, "TcgaTargetGtex_gene_expected_count.gz")
PHENO_FILE <- file.path(BASE_DIR, "TcgaTargetGTEX_phenotype.txt.gz")

dir.create(RES_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(META_DIR, showWarnings = FALSE, recursive = TRUE)

WORKERS     <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8"))
CODING_ONLY <- TRUE

register(MulticoreParam(workers = WORKERS, progressbar = FALSE))
setDTthreads(WORKERS)

cat("=================================================================\n")
cat(" Script 01: DESeq2 — Glioma vs Healthy Normal Brain (Xena/Toil)\n")
cat(sprintf(" Base directory : %s\n", BASE_DIR))
cat(sprintf(" Workers        : %d\n", WORKERS))
cat(sprintf(" Coding only    : %s\n", CODING_ONLY))
cat("=================================================================\n\n")

for (f in c(COUNT_FILE, PHENO_FILE))
  if (!file.exists(f)) stop(sprintf("Input file not found: %s", f))

# =============================================================================
# STEP 1: Load phenotype and select samples
# =============================================================================
cat("--- Step 1: Loading phenotype annotations ---\n")

pheno <- fread(PHENO_FILE, sep = "\t", header = TRUE, data.table = TRUE)
setnames(pheno, gsub("^_", "",  names(pheno)))
setnames(pheno, gsub(" ",  ".", names(pheno)))
cat(sprintf("  Phenotype table: %d samples\n", nrow(pheno)))

# Diagnostic: print ALL GTEx Brain tissue labels found in the file
# so we can confirm our three target regions are present
gtex_brain_labels <- sort(unique(
  pheno[study == "GTEX" & grepl("Brain", primary.disease.or.tissue,
                                 ignore.case = TRUE),
        primary.disease.or.tissue]))
cat("  All GTEx Brain tissue labels in file:\n")
for (l in gtex_brain_labels) cat(sprintf("    %s\n", l))
cat("\n")

# ---------------------------------------------------------------------------
# GTEx: THREE brain regions
#   Frontal Cortex (BA9)               <- primary glioma-matched region
#   Anterior Cingulate Cortex (BA24)   <- adjacent cortical region
#   Cortex (broad)                     <- general cortex
#
# These exactly match the GCT files:
#   gene_reads_v11_brain_frontal_cortex_ba9.gct
#   gene_reads_v11_brain_anterior_cingulate_cortex_ba24.gct
#   gene_reads_v11_brain_cortex.gct
# ---------------------------------------------------------------------------
gtex_ba9  <- pheno[
  study == "GTEX" &
  grepl("Frontal Cortex", primary.disease.or.tissue, ignore.case = TRUE) &
  grepl("BA9",            primary.disease.or.tissue, ignore.case = TRUE),
  sample]

gtex_ba24 <- pheno[
  study == "GTEX" &
  grepl("Anterior Cingulate", primary.disease.or.tissue, ignore.case = TRUE) &
  grepl("BA24",               primary.disease.or.tissue, ignore.case = TRUE),
  sample]

gtex_ctx  <- pheno[
  study == "GTEX" &
  grepl("^Brain - Cortex$", primary.disease.or.tissue, ignore.case = TRUE),
  sample]

gtex_keep <- unique(c(gtex_ba9, gtex_ba24, gtex_ctx))

cat(sprintf("  GTEx Brain - Frontal Cortex (BA9)              : %d samples\n",
            length(gtex_ba9)))
cat(sprintf("  GTEx Brain - Anterior Cingulate Cortex (BA24)  : %d samples\n",
            length(gtex_ba24)))
cat(sprintf("  GTEx Brain - Cortex                            : %d samples\n",
            length(gtex_ctx)))
cat(sprintf("  GTEx total (3 regions combined)                : %d samples\n\n",
            length(gtex_keep)))

# Warn if any region returned zero samples
if (length(gtex_ba9)  == 0L) warning("GTEx BA9 returned 0 samples — check label above.")
if (length(gtex_ba24) == 0L) warning("GTEx BA24 returned 0 samples — check label above.")
if (length(gtex_ctx)  == 0L) warning("GTEx Cortex returned 0 samples — check label above.")
if (length(gtex_keep) == 0L) stop("No GTEx brain samples found at all.")

# ---------------------------------------------------------------------------
# TCGA: GBM + LGG, primary tumors only
# Diagnostic: print exact values so we know what to match against
# ---------------------------------------------------------------------------
cat("  Unique 'study' column values:\n")
for (s in sort(unique(pheno$study))) cat(sprintf("    '%s'\n", s))
cat("\n")

cat("  Unique 'sample_type' values:\n")
for (s in sort(unique(pheno$sample_type))) cat(sprintf("    '%s'\n", s))
cat("\n")

cat("  All TCGA primary.disease.or.tissue labels:\n")
for (l in sort(unique(pheno[study == "TCGA", primary.disease.or.tissue])))
  cat(sprintf("    '%s'\n", l))
cat("\n")

cat("  All TCGA primary_site values:\n")
if ("primary_site" %in% names(pheno)) {
  for (l in sort(unique(pheno[study == "TCGA", primary_site])))
    cat(sprintf("    '%s'\n", l))
} else {
  cat("    (column not present)\n")
}
cat("\n")

# Filter strategy: use primary_site == "Brain" to catch both GBM and LGG
# regardless of how the disease name is spelled in the file.
# Falls back to disease name matching if primary_site is not available.
if ("primary_site" %in% names(pheno)) {
  tcga_keep <- pheno[
    study == "TCGA" &
    primary_site == "Brain" &
    grepl("Primary Tumor|Primary$", sample_type, ignore.case = TRUE),
    sample]
  cat(sprintf("  Filter method: primary_site == 'Brain'\n"))
} else {
  # Fallback: match disease name
  tcga_keep <- pheno[
    study == "TCGA" &
    grepl("glioblastoma|glioma|GBM|LGG|lower.grade", primary.disease.or.tissue,
          ignore.case = TRUE) &
    grepl("Primary Tumor|Primary$", sample_type, ignore.case = TRUE),
    sample]
  cat(sprintf("  Filter method: disease name grep (primary_site column absent)\n"))
}

# Show breakdown by disease label for transparency
if (length(tcga_keep) > 0L) {
  sub_labels <- unique(pheno[sample %in% tcga_keep, primary.disease.or.tissue])
  cat("  Disease labels included:\n")
  for (l in sort(sub_labels)) cat(sprintf("    '%s'\n", l))
  cat("\n")
}

cat(sprintf("  TCGA GBM + LGG (Primary Tumor) : %d samples\n", length(tcga_keep)))
cat(sprintf("  Grand total                     : %d samples\n\n",
            length(gtex_keep) + length(tcga_keep)))

if (length(tcga_keep) == 0L)
  stop("No TCGA GBM/LGG Primary Tumor samples found.")

# =============================================================================
# STEP 2: Load count matrix — selected columns only
#
# WORKAROUND FOR R's 2^31-1 BYTE STRING LIMIT
#
# The full Xena matrix has ~20,000 tab-separated columns per row.
# When fread decompresses the .gz file, each row becomes a single character
# string that exceeds R's 2 GB limit, even with select = <indices>.
#
# Solution: use zcat + awk on the command line to extract ONLY the ~950
# columns we need into a small temporary TSV, then fread that file.
# =============================================================================
cat("--- Step 2: Loading count matrix (selected columns only) ---\n")
cat("  Reading header...\n")

con      <- gzcon(file(COUNT_FILE, "rb"))
hdr_line <- readLines(con, n = 1L)
close(con)
all_cols <- strsplit(hdr_line, "\t")[[1]]

all_want    <- c(gtex_keep, tcga_keep)
matched_idx <- which(all_cols %in% all_want)
missing     <- setdiff(all_want, all_cols[matched_idx])
if (length(missing) > 0L)
  warning(sprintf("%d sample IDs not found in count matrix (first 5: %s)",
                  length(missing), paste(head(missing, 5L), collapse = ", ")))

cat(sprintf("  Loading %d matched columns...\n", length(matched_idx)))

# Build awk column extraction command
# Column 1 = gene IDs, then only the matched sample columns
col_indices <- c(1L, matched_idx)

# Write an awk script to a file — avoids shell command-line length limits
# when there are ~950 columns to extract
AWK_SCRIPT <- file.path(BASE_DIR, "tmp_extract_cols.awk")
awk_body <- paste0("$", col_indices, collapse = ",")
writeLines(
  c(
    "BEGIN { FS=\"\\t\"; OFS=\"\\t\" }",
    sprintf("{ print %s }", awk_body)
  ),
  AWK_SCRIPT
)

SUBSET_FILE <- file.path(BASE_DIR, "tmp_count_subset.tsv")
awk_cmd <- sprintf(
  "zcat '%s' | awk -f '%s' > '%s'",
  COUNT_FILE, AWK_SCRIPT, SUBSET_FILE
)

cat("  Extracting columns via zcat | awk (this may take 10–20 min)...\n")
t0  <- proc.time()
# Use bash -c so pipefail catches real errors; ignore awk SIGPIPE (141)
ret <- system(sprintf("bash -c 'set -o pipefail; %s'", awk_cmd))
if (ret != 0L && ret != 141L) {
  # Verify the output file exists and is non-empty as a fallback check
  if (!file.exists(SUBSET_FILE) || file.size(SUBSET_FILE) < 1000L)
    stop("awk column extraction failed — check disk space and file integrity.")
}
cat(sprintf("  awk extraction done: %.1f min\n", (proc.time() - t0)[["elapsed"]] / 60))

# Verify the temp file looks reasonable
subset_size <- file.size(SUBSET_FILE)
cat(sprintf("  Subset file size: %.1f MB\n", subset_size / 1e6))
if (subset_size < 1e6)
  stop("Subset file suspiciously small — extraction likely failed.")

t0     <- proc.time()
raw_dt <- fread(SUBSET_FILE, sep = "\t", header = TRUE, data.table = TRUE)
cat(sprintf("  Loaded in %.1f min: %d genes x %d samples\n",
            (proc.time() - t0)[["elapsed"]] / 60,
            nrow(raw_dt), ncol(raw_dt) - 1L))

# Clean up temp files
unlink(SUBSET_FILE)
unlink(AWK_SCRIPT)

# The first column is the gene ID column — rename it defensively
# (Xena calls it "sample" but that could change)
setnames(raw_dt, 1L, "gene_id")

# Strip Ensembl version suffix and deduplicate
raw_dt[, gene_id := sub("\\.[0-9]+$", "", gene_id)]
raw_dt <- raw_dt[!duplicated(gene_id)]
gene_ids <- raw_dt$gene_id
raw_dt[, gene_id := NULL]

# Validate all remaining columns are numeric before back-transform
col_classes <- sapply(raw_dt, class)
non_numeric <- names(col_classes)[!col_classes %in% c("numeric", "integer", "double")]
if (length(non_numeric) > 0L) {
  cat(sprintf("  WARNING: %d non-numeric columns detected — coercing.\n",
              length(non_numeric)))
  cat(sprintf("  First 5: %s\n", paste(head(non_numeric, 5), collapse = ", ")))
  for (col in non_numeric)
    set(raw_dt, j = col, value = as.numeric(raw_dt[[col]]))
}

# Back-transform log2(x+1) -> integer counts
cat("  Back-transforming log2(x+1) -> round(2^x - 1)...\n")
mat <- as.matrix(raw_dt)
rm(raw_dt); gc()

# Check for NAs introduced by coercion
n_na <- sum(is.na(mat))
if (n_na > 0L) {
  cat(sprintf("  WARNING: %d NAs in matrix — setting to 0.\n", n_na))
  mat[is.na(mat)] <- 0
}

count_mat <- round(2^mat - 1)
rm(mat); gc()
count_mat[count_mat < 0L] <- 0L
storage.mode(count_mat) <- "integer"
rownames(count_mat) <- gene_ids

cat(sprintf("  Count matrix: %d genes x %d samples\n\n",
            nrow(count_mat), ncol(count_mat)))

# =============================================================================
# STEP 3: ENSG -> symbol map
# =============================================================================
cat("--- Step 3: Gene annotation (ENSG -> symbol) ---\n")

ensg_map_path <- file.path(META_DIR, "ensg_to_symbol.rds")

if (file.exists(ensg_map_path)) {

  ensg_map <- readRDS(ensg_map_path)
  cat(sprintf("  Loaded cached map: %d entries\n", nrow(ensg_map)))

} else if (requireNamespace("biomaRt", quietly = TRUE)) {

  cat("  Fetching from Ensembl biomaRt (requires internet)...\n")
  suppressPackageStartupMessages(library(biomaRt))
  mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",
                     version = 100)
  bm   <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol",
                                "gene_biotype"), mart = mart)
  ensg_map <- as.data.table(bm)
  setnames(ensg_map, c("gene_id", "symbol", "gene_type"))
  ensg_map <- ensg_map[symbol != ""]
  ensg_map <- unique(ensg_map, by = "gene_id")
  saveRDS(ensg_map, ensg_map_path)
  cat(sprintf("  Fetched and cached: %d entries\n", nrow(ensg_map)))

} else {

  cat("  WARNING: biomaRt unavailable — using ENSG IDs as symbols.\n")
  cat("  Install biomaRt and delete the cache to get real symbols.\n")
  ensg_map <- data.table(gene_id   = rownames(count_mat),
                         symbol    = rownames(count_mat),
                         gene_type = "unknown")
  saveRDS(ensg_map, ensg_map_path)
}

# =============================================================================
# STEP 4: Protein-coding filter
# =============================================================================
cat("\n--- Step 4: Gene filtering ---\n")
cat(sprintf("  Genes before filter: %d\n", nrow(count_mat)))

if (CODING_ONLY && "gene_type" %in% names(ensg_map) &&
    any(ensg_map$gene_type == "protein_coding")) {
  coding_ids <- intersect(ensg_map[gene_type == "protein_coding", gene_id],
                          rownames(count_mat))
  count_mat  <- count_mat[coding_ids, , drop = FALSE]
  cat(sprintf("  After protein-coding filter: %d genes\n", nrow(count_mat)))
  nat8l_id <- ensg_map[symbol == "NAT8L", gene_id][1L]
  cat(sprintf("  Sanity check — NAT8L present: %s\n",
              ifelse(!is.na(nat8l_id) && nat8l_id %in% rownames(count_mat),
                     "YES [OK]", "not found")))
} else {
  cat("  Skipping (CODING_ONLY=FALSE or gene_type unavailable)\n")
}

# =============================================================================
# STEP 5: Build sample metadata
# =============================================================================
cat("\n--- Step 5: Building sample metadata ---\n")

final_samples <- colnames(count_mat)
tcga_final    <- intersect(final_samples, tcga_keep)
gtex_final    <- intersect(final_samples, gtex_keep)

# Assign GTEx sub-region label (for logging; not used in model)
gtex_region <- fcase(
  gtex_final %in% gtex_ba9,  "BA9",
  gtex_final %in% gtex_ba24, "BA24",
  gtex_final %in% gtex_ctx,  "Cortex",
  default = "Brain_other"
)

# Pull GBM vs LGG label from phenotype file
pheno_sub <- pheno[sample %in% tcga_final,
                   .(sample, primary.disease.or.tissue)]
pheno_sub[, cohort := ifelse(grepl("GBM", primary.disease.or.tissue,
                                    ignore.case = TRUE), "GBM", "LGG")]

metadata <- rbind(
  data.table(sample_id = pheno_sub$sample,
             condition = "glioma",
             cohort    = pheno_sub$cohort),
  data.table(sample_id = gtex_final,
             condition = "normal",
             cohort    = paste0("GTEx_", gtex_region))
)
metadata[, condition := factor(condition, levels = c("normal", "glioma"))]

cat(sprintf("  Total samples  : %d\n", nrow(metadata)))
cat(sprintf("  Glioma (TCGA)  : %d  [GBM: %d | LGG: %d]\n",
            sum(metadata$condition == "glioma"),
            sum(metadata$cohort == "GBM"),
            sum(metadata$cohort == "LGG")))
cat(sprintf("  Normal (GTEx)  : %d  [BA9: %d | BA24: %d | Cortex: %d]\n",
            sum(metadata$condition == "normal"),
            sum(metadata$cohort == "GTEx_BA9"),
            sum(metadata$cohort == "GTEx_BA24"),
            sum(metadata$cohort == "GTEx_Cortex")))

count_mat <- count_mat[, metadata$sample_id, drop = FALSE]
stopifnot(identical(colnames(count_mat), metadata$sample_id))

# =============================================================================
# STEP 6: Low-count filtering
# =============================================================================
cat("\n--- Step 6: Low-count filtering ---\n")

min_grp   <- min(table(metadata$condition))
keep      <- rowSums(count_mat >= 10L) >= min_grp
cat(sprintf("  Min group size: %d  |  Genes: %d -> %d\n",
            min_grp, nrow(count_mat), sum(keep)))
count_mat <- count_mat[keep, , drop = FALSE]

# =============================================================================
# STEP 7: DESeq2
# =============================================================================
cat("\n--- Step 7: DESeq2 (glioma vs normal) ---\n")
cat(sprintf("  %d samples x %d genes\n", ncol(count_mat), nrow(count_mat)))

meta_df           <- as.data.frame(metadata)
rownames(meta_df) <- meta_df$sample_id

dds           <- DESeqDataSetFromMatrix(count_mat, meta_df, ~ condition)
dds$condition <- relevel(dds$condition, ref = "normal")

cat("  Running DESeq2...\n")
t0  <- proc.time()
dds <- DESeq(dds, parallel = TRUE, BPPARAM = MulticoreParam(WORKERS))
cat(sprintf("  Done: %.1f min\n", (proc.time() - t0)[["elapsed"]] / 60))

res    <- results(dds,
                  contrast = c("condition", "glioma", "normal"),
                  alpha    = 0.05,
                  parallel = TRUE,
                  BPPARAM  = MulticoreParam(WORKERS))
res_dt <- as.data.table(as.data.frame(res), keep.rownames = "gene_id")
res_dt <- merge(res_dt, ensg_map[, .(gene_id, symbol)],
                by = "gene_id", all.x = TRUE)

n_sig <- sum(res_dt$padj < 0.05 & !is.na(res_dt$padj))
cat(sprintf("  DEGs (padj<0.05): %d  [Up: %d | Down: %d]\n",
            n_sig,
            sum(res_dt$padj < 0.05 & res_dt$log2FoldChange > 0, na.rm = TRUE),
            sum(res_dt$padj < 0.05 & res_dt$log2FoldChange < 0, na.rm = TRUE)))

# =============================================================================
# STEP 8: Save outputs
# =============================================================================
cat("\n--- Step 8: Saving outputs ---\n")

fwrite(res_dt,
       file.path(RES_DIR, "deseq2_results_glioma_vs_normal.tsv"),
       sep = "\t")
saveRDS(dds,
        file.path(RES_DIR, "DESeq2_dds_glioma_vs_normal.rds"))

# Ranked gene list for fgsea (Wald stat: >0 = higher in glioma)
ranked <- sort(setNames(res_dt[!is.na(stat), stat],
                         res_dt[!is.na(stat), gene_id]),
               decreasing = TRUE)

# NAT8L sanity check (expected depleted in glioma vs normal brain)
nat8l_id <- ensg_map[symbol == "NAT8L", gene_id][1L]
if (!is.na(nat8l_id) && nat8l_id %in% names(ranked))
  cat(sprintf("  NAT8L stat: %.3f  %s\n", ranked[nat8l_id],
              ifelse(ranked[nat8l_id] < 0,
                     "[EXPECTED: depleted in glioma]",
                     "[WARNING: check direction]")))

saveRDS(ranked, file.path(RES_DIR, "glioma_vs_normal_ranked_genes.rds"))

cat(sprintf("  deseq2_results_glioma_vs_normal.tsv\n"))
cat(sprintf("  DESeq2_dds_glioma_vs_normal.rds\n"))
cat(sprintf("  glioma_vs_normal_ranked_genes.rds\n"))
cat(sprintf("  Ranked list: %d genes | range [%.2f, %.2f]\n",
            length(ranked), min(ranked), max(ranked)))

# Session info
sink(file.path(RES_DIR, "01_deseq2_session_info.txt"))
cat("DESeq2 Glioma vs Healthy Brain — Session Info\n")
cat(format(Sys.time(), "Run: %Y-%m-%d %H:%M:%S %Z\n\n"))
cat(sprintf("Counts  : %s\n", COUNT_FILE))
cat(sprintf("Pheno   : %s\n", PHENO_FILE))
cat(sprintf("GTEx regions: BA9=%d | BA24=%d | Cortex=%d | Total=%d\n",
            length(gtex_ba9), length(gtex_ba24), length(gtex_ctx),
            length(gtex_keep)))
cat(sprintf("TCGA    : %d samples (GBM + LGG, Primary Tumor)\n",
            length(tcga_keep)))
cat(sprintf("Genes   : %d (protein-coding: %s)\n",
            nrow(count_mat), CODING_ONLY))
cat(sprintf("Workers : %d\n\n", WORKERS))
print(sessionInfo())
sink()

cat("\n=================================================================\n")
cat(" Script 01 complete\n")
cat(sprintf(" Results: %s\n", RES_DIR))
cat("=================================================================\n")# Filter to samples with IDH annotation
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
