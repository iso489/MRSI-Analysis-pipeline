#!/usr/bin/env Rscript
# =============================================================================
# 05 - Donor-level reanalysis: the PRIMARY transcriptomic inference
# =============================================================================
#
# WHY THIS EXISTS
# The published comparison treats 289 GTEx cortical SAMPLES as 289 independent
# normal brains. They are not. The 289 samples come from 164 unique GTEx
# donors: 73 donors contribute one sample, 59 contribute two, 31 contribute
# three and one (GTEX-NPJ8) contributes five. 91 donors therefore contribute
# repeated observations, because GTEx sampled up to three cortical regions
# (Frontal Cortex BA9, Anterior Cingulate BA24, Cortex) from the same brain.
#
# The published model is ~condition with no donor term, and the GSVA
# comparison is an ordinary Wilcoxon rank-sum test. Neither accounts for the
# repeated sampling, so the nominal P values, the BH-adjusted values and the
# DESeq2 Wald statistics that feed the pre-ranked GSEA are all anti-conservative
# to an unknown degree. The size-matched random-gene-set calibration does not
# repair this: it re-uses the same dependent ranking.
#
# WHAT THIS SCRIPT DOES
#   0. Emits the exact sample -> donor -> region manifest (publishable).
#   1. A0  reproduces the published sample-level analysis (regression check).
#   2. A1  PRIMARY: within-donor aggregation. GTEx counts are summed across a
#          donor's eligible cortical samples (summing is the correct pooling for
#          negative-binomial counts), giving 164 independent normal donors vs
#          662 TCGA participants (one primary tumor each - verified, no repeats).
#   3. A2  SENSITIVITY: one prespecified sample per donor, region priority
#          BA9 > BA24 > Cortex, ties broken by the lexicographically smallest
#          sample ID. Deterministic, no random draw.
#   4. A3  SENSITIVITY: region-specific analyses, each within-donor aggregated
#          (BA9 99 donors, BA24 81 donors, Cortex 102 donors).
#   5. Re-runs fgsea, GSVA + Wilcoxon, NAT8L and the 500-set size-matched
#      random calibration on the primary donor-level fit.
#
# WHAT IT LICENSES
#   LICENSES: donor-level inferential quantities for the external transcriptomic
#             context analysis.
#   DOES NOT: remove the study-source confounding (surgical tumor vs postmortem
#             cortex), make the comparison spatial, patient-matched or causal.
#
# Run:  Rscript 05_donor_level_transcriptomics.R
#       (set GVN_BASE to override /scratch/$USER/glioma_vs_normal)
# =============================================================================

suppressMessages({
  library(DESeq2); library(fgsea); library(GSVA); library(data.table)
})

USER <- Sys.getenv("USER")
  grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
BASE_DIR <- Sys.getenv("GVN_BASE", unset = file.path("/scratch", USER, "glioma_vs_normal"))

DDS_RDS  <- file.path(BASE_DIR, "results", "DESeq2_dds_glioma_vs_normal.rds")
GMT      <- file.path(BASE_DIR, "metadata", "pathways.gmt")
SYM_RDS  <- file.path(BASE_DIR, "metadata", "ensg_to_symbol.rds")
PC_RES   <- file.path(BASE_DIR, "results", "deseq2_results_protein_coding.tsv")
OUT      <- file.path(BASE_DIR, "results", "donor_level")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

for (p in c(DDS_RDS, GMT, SYM_RDS, PC_RES))
  if (!file.exists(p)) stop("Input not found: ", p)

SEED <- 20260814
set.seed(SEED)
B_RANDOM <- 500L

log <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

# =============================================================================
# 0. Load, restrict to protein-coding, build the donor manifest
# =============================================================================
log("loading dds")
dds0 <- readRDS(DDS_RDS)
log("dds: ", nrow(dds0), " genes x ", ncol(dds0), " samples")

pc <- fread(PC_RES, select = "gene_id")
keep <- rownames(dds0) %in% pc$gene_id
log("protein-coding retained: ", sum(keep), " / ", nrow(dds0))
cat("protein-coding genes:", sum(keep), "\n")

cnt  <- counts(dds0)[keep, , drop = FALSE]
cd   <- as.data.frame(colData(dds0))[, c("sample_id", "condition", "cohort")]
cd$condition <- as.character(cd$condition)
cd$cohort    <- as.character(cd$cohort)

# Donor key: GTEx sample IDs are GTEX-<donor>-<tissue site>-SM-<aliquot>, so the
# first two hyphen-separated fields identify the brain. TCGA IDs are
# TCGA-<TSS>-<participant>-<sample>, so the first three identify the patient.
donor_of <- function(id) {
  n <- ifelse(startsWith(id, "GTEX"), 2L, 3L)
  vapply(seq_along(id), function(i)
    paste(strsplit(id[i], "-", fixed = TRUE)[[1]][seq_len(n[i])], collapse = "-"),
    character(1))
}
cd$donor  <- donor_of(cd$sample_id)
cd$region <- ifelse(cd$condition == "normal", sub("^GTEx_", "", cd$cohort), NA_character_)

gt <- cd[cd$condition == "normal", ]
tc <- cd[cd$condition == "glioma", ]
log("GTEx: ", nrow(gt), " samples from ", length(unique(gt$donor)), " donors")
log("TCGA: ", nrow(tc), " samples from ", length(unique(tc$donor)), " participants")
stopifnot(nrow(gt) == 289L, length(unique(gt$donor)) == 164L,
          nrow(tc) == 662L, length(unique(tc$donor)) == 662L)

spd <- table(table(gt$donor))
log("GTEx samples per donor: ",
    paste(sprintf("%s sample(s): %d donor(s)", names(spd), as.integer(spd)), collapse = "; "))

mf <- data.table(cd)
mf[, n_samples_from_donor := .N, by = donor]
setcolorder(mf, c("sample_id", "donor", "condition", "region", "n_samples_from_donor"))
fwrite(mf[order(condition, donor, sample_id)],
       file.path(OUT, "sample_donor_manifest.tsv"), sep = "\t")
log("wrote sample_donor_manifest.tsv")

# =============================================================================
# helpers
# =============================================================================
run_de <- function(counts_mat, col_data, label) {
  d <- DESeqDataSetFromMatrix(countData = counts_mat,
                              colData   = col_data,
                              design    = ~condition)
  d$condition <- relevel(factor(d$condition), ref = "normal")
  d <- DESeq(d, quiet = TRUE)
  r <- results(d, contrast = c("condition", "glioma", "normal"))
  n_de <- sum(r$padj < 0.05, na.rm = TRUE); n_t <- sum(!is.na(r$padj))
  log(label, ": n=", ncol(d), " (", sum(col_data$condition == "normal"), " normal / ",
      sum(col_data$condition == "glioma"), " glioma)  DE=", n_de, "/", n_t,
      " = ", round(100 * n_de / n_t, 1), "%")
  list(dds = d, res = r)
}

sym <- readRDS(SYM_RDS)
ranked <- function(res) {
  dt <- as.data.table(as.data.frame(res), keep.rownames = "gene_id")
  dt[, symbol := sym$symbol[match(gene_id, sym$gene_id)]]
  rk <- dt[!is.na(stat) & !is.na(symbol) & symbol != ""]
  rk <- rk[order(-abs(rk$stat))][!duplicated(symbol)]
  list(dt = dt, ranks = sort(setNames(rk$stat, rk$symbol), decreasing = TRUE))
}

pw <- gmtPathways(GMT)

run_gsea <- function(ranks, label) {
  fg <- fgsea(pathways = pw, stats = ranks, minSize = 5, maxSize = 500,
              nPermSimple = 10000)
  fg <- fg[order(padj)]
  fg[, analysis := label]
  fg[, .(analysis, pathway, NES, pval, padj, size)]
}

run_gsva <- function(dds_fit, label) {
  vsd <- vst(dds_fit, blind = FALSE)
  m   <- assay(vsd)
  rownames(m) <- sym$symbol[match(rownames(m), sym$gene_id)]
  m   <- m[!is.na(rownames(m)) & rownames(m) != "", ]
  agg <- rowsum(m, rownames(m))
  m   <- agg / as.vector(table(rownames(m))[rownames(agg)])
  gs  <- gsva(gsvaParam(as.matrix(m), pw, kcdf = "Gaussian"), verbose = FALSE)
  grp <- as.character(colData(dds_fit)$condition)
  st  <- rbindlist(lapply(rownames(gs), function(p) {
    a <- gs[p, grp == "glioma"]; b <- gs[p, grp == "normal"]
    data.table(analysis = label, pathway = p,
               mean_glioma = mean(a), mean_normal = mean(b),
               delta = mean(a) - mean(b),
               p = wilcox.test(a, b)$p.value)
  }))
  st[, padj := p.adjust(p, method = "BH")]
  list(stats = st[order(padj)], scores = gs)
}

genes_of_interest <- c("NAT8L", "SNAP25", "NEFL", "SYN1", "SLC12A5", "GFAP")
single_genes <- function(dt, label) {
  out <- dt[symbol %in% genes_of_interest,
            .(analysis = label, symbol, log2FoldChange, stat, pvalue, padj)]
  out[order(match(symbol, genes_of_interest))]
}

# =============================================================================
# 1. A0 - reproduce the published SAMPLE-level analysis (regression check)
# =============================================================================
log("=== A0: published sample-level analysis (289 samples vs 662) ===")
a0 <- run_de(cnt, cd[, c("sample_id", "condition", "cohort")], "A0 sample-level")
r0 <- ranked(a0$res)

# =============================================================================
# 2. A1 - PRIMARY donor-level analysis (within-donor count aggregation)
# =============================================================================
log("=== A1 PRIMARY: within-donor aggregation (164 donors vs 662) ===")
cnt_a1 <- t(rowsum(t(cnt), cd$donor))                 # genes x donors
cd_a1  <- unique(data.table(cd)[, .(donor, condition)])
setkey(cd_a1, donor)
cd_a1  <- as.data.frame(cd_a1[colnames(cnt_a1)])
rownames(cd_a1) <- cd_a1$donor
cd_a1$sample_id <- cd_a1$donor
stopifnot(ncol(cnt_a1) == 826L, sum(cd_a1$condition == "normal") == 164L)

a1 <- run_de(cnt_a1, cd_a1, "A1 donor-level (PRIMARY)")
r1 <- ranked(a1$res)
fwrite(r1$dt, file.path(OUT, "A1_deseq2_donor_level.tsv"), sep = "\t")

# =============================================================================
# 3. A2 - one prespecified sample per donor
# =============================================================================
log("=== A2 SENSITIVITY: one sample per donor (BA9 > BA24 > Cortex) ===")
prio <- c(BA9 = 1L, BA24 = 2L, Cortex = 3L)
gt_dt <- data.table(gt)
gt_dt[, rank := prio[region]]
pick <- gt_dt[order(donor, rank, sample_id), .SD[1L], by = donor]$sample_id
stopifnot(length(pick) == 164L)
sel_a2 <- c(pick, tc$sample_id)
a2 <- run_de(cnt[, sel_a2, drop = FALSE],
             cd[sel_a2, c("sample_id", "condition", "cohort")],
             "A2 one-per-donor")
r2 <- ranked(a2$res)
fwrite(data.table(sample_id = pick), file.path(OUT, "A2_selected_samples.tsv"), sep = "\t")

# =============================================================================
# 4. A3 - region-specific, within-donor aggregated
# =============================================================================
log("=== A3 SENSITIVITY: region-specific ===")
a3_res <- list()
for (rg in c("BA9", "BA24", "Cortex")) {
  ids <- gt$sample_id[gt$region == rg]
  sub_cd  <- cd[c(ids, tc$sample_id), ]
  sub_cnt <- cnt[, c(ids, tc$sample_id), drop = FALSE]
  m  <- t(rowsum(t(sub_cnt), sub_cd$donor))
  c2 <- unique(data.table(sub_cd)[, .(donor, condition)]); setkey(c2, donor)
  c2 <- as.data.frame(c2[colnames(m)]); rownames(c2) <- c2$donor
  c2$sample_id <- c2$donor
  fit <- run_de(m, c2, paste0("A3 ", rg))
  a3_res[[rg]] <- ranked(fit$res)
}

# =============================================================================
# 5. Downstream on the PRIMARY fit: fgsea, GSVA, single genes, calibration
# =============================================================================
log("=== enrichment across analyses ===")
gsea_all <- rbindlist(list(
  run_gsea(r0$ranks, "A0 sample-level"),
  run_gsea(r1$ranks, "A1 donor-level (PRIMARY)"),
  run_gsea(r2$ranks, "A2 one-per-donor"),
  rbindlist(lapply(names(a3_res), function(rg) run_gsea(a3_res[[rg]]$ranks, paste0("A3 ", rg))))
))
fwrite(gsea_all, file.path(OUT, "fgsea_all_analyses.tsv"), sep = "\t")

log("GSVA")
g0 <- run_gsva(a0$dds, "A0 sample-level")
g1 <- run_gsva(a1$dds, "A1 donor-level (PRIMARY)")
gsva_all <- rbindlist(list(g0$stats, g1$stats))
fwrite(gsva_all, file.path(OUT, "gsva_all_analyses.tsv"), sep = "\t")
write.csv(g1$scores, file.path(OUT, "A1_gsva_scores_donor_level.csv"))

sg <- rbindlist(list(
  single_genes(r0$dt, "A0 sample-level"),
  single_genes(r1$dt, "A1 donor-level (PRIMARY)"),
  single_genes(r2$dt, "A2 one-per-donor"),
  rbindlist(lapply(names(a3_res), function(rg) single_genes(a3_res[[rg]]$dt, paste0("A3 ", rg))))
))
fwrite(sg, file.path(OUT, "single_genes_all_analyses.tsv"), sep = "\t")

log("size-matched random-set calibration on the primary donor-level ranking")
ranks1 <- r1$ranks
obs <- fgsea(pw, ranks1, minSize = 5, maxSize = 500, nPermSimple = 10000)[order(padj)]
universe <- names(ranks1)
calib <- rbindlist(lapply(seq_len(nrow(obs)), function(i) {
  nm <- obs$pathway[i]; sz <- obs$size[i]
  rnd <- setNames(lapply(seq_len(B_RANDOM), function(b) sample(universe, sz)),
                  sprintf("RND_%s_%04d", sz, seq_len(B_RANDOM)))
  fr <- suppressWarnings(fgsea(rnd, ranks1, minSize = 5, maxSize = 500,
                               nPermSimple = 1000))
  data.table(pathway = nm, size = sz,
             observed_NES = obs$NES[i], observed_padj = obs$padj[i],
             random_NES_q025 = as.numeric(quantile(fr$NES, .025, na.rm = TRUE)),
             random_NES_q975 = as.numeric(quantile(fr$NES, .975, na.rm = TRUE)),
             random_absNES_max = max(abs(fr$NES), na.rm = TRUE),
             n_random = nrow(fr),
             calibrated_p = (sum(abs(fr$NES) >= abs(obs$NES[i]), na.rm = TRUE) + 1) /
                            (nrow(fr) + 1),
             random_frac_FDR05 = mean(fr$padj < 0.05, na.rm = TRUE))
}))
calib <- calib[order(calibrated_p, -abs(observed_NES))]
fwrite(calib, file.path(OUT, "A1_negative_control_gsea.tsv"), sep = "\t")

# =============================================================================
# 6. Summary
# =============================================================================
cat("\n================ DONOR-LEVEL SUMMARY ================\n")
cat("\n-- Differential expression --\n")
de_summary <- rbindlist(lapply(
  list(list("A0 sample-level", a0$res, 289L),
       list("A1 donor-level (PRIMARY)", a1$res, 164L),
       list("A2 one-per-donor", a2$res, 164L)),
  function(x) data.table(analysis = x[[1]], n_normal = x[[3]],
    n_DE = sum(x[[2]]$padj < 0.05, na.rm = TRUE),
    n_tested = sum(!is.na(x[[2]]$padj)),
    pct = round(100 * sum(x[[2]]$padj < 0.05, na.rm = TRUE) / sum(!is.na(x[[2]]$padj)), 1))))
print(de_summary)
fwrite(de_summary, file.path(OUT, "de_summary.tsv"), sep = "\t")

cat("\n-- fgsea, primary vs published --\n")
cmp <- merge(gsea_all[analysis == "A0 sample-level", .(pathway, NES_A0 = NES, padj_A0 = padj)],
             gsea_all[analysis == "A1 donor-level (PRIMARY)", .(pathway, NES_A1 = NES, padj_A1 = padj)],
             by = "pathway")
print(cmp[order(NES_A1), .(pathway, NES_A0 = round(NES_A0, 2), NES_A1 = round(NES_A1, 2),
                           padj_A0 = signif(padj_A0, 3), padj_A1 = signif(padj_A1, 3))])

cat("\n-- GSVA, primary vs published --\n")
gcmp <- merge(gsva_all[analysis == "A0 sample-level", .(pathway, d_A0 = delta, p_A0 = padj)],
              gsva_all[analysis == "A1 donor-level (PRIMARY)", .(pathway, d_A1 = delta, p_A1 = padj)],
              by = "pathway")
print(gcmp[order(d_A1), .(pathway, d_A0 = round(d_A0, 3), d_A1 = round(d_A1, 3),
                          padj_A0 = signif(p_A0, 3), padj_A1 = signif(p_A1, 3))])

cat("\n-- NAT8L across analyses --\n")
print(sg[symbol == "NAT8L", .(analysis, log2FoldChange = round(log2FoldChange, 3),
                              stat = round(stat, 1), padj = signif(padj, 3))])

cat("\n-- calibration (primary) --\n")
print(calib[, .(pathway, size, NES = round(observed_NES, 2),
                calibrated_p = signif(calibrated_p, 3))])

writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo.txt"))
log("done -> ", OUT)
