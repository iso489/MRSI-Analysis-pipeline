#!/usr/bin/env Rscript
# =============================================================================
# 02 - Size-matched random gene sets as negative controls for the enrichment
#      analysis
# =============================================================================
#
# WHY THIS ANALYSIS
# 92.1% of tested protein-coding genes are differentially expressed between
# bulk glioma and normal cortex. Both reviewers raised the obvious objection:
# with a genome-wide shift that large, any gene set will look enriched, so the
# prespecified MRS-relevant results may be uninformative. This tests that
# objection empirically rather than arguing about it.
#
# DESIGN
# For each prespecified gene set we draw B random gene sets of EXACTLY the same
# size from the same tested universe (the ranked gene list actually used) and
# compute enrichment for each. The reported calibration is the empirical
# proportion of size-matched random sets reaching an |NES| at least as extreme
# as the observed one - a competitive gene-set-level null.
#
# WHY NOT SIMPLY ADD THE RANDOM SETS TO THE MAIN fgsea CALL
# That would inflate the Benjamini-Hochberg denominator and artificially
# depress the prespecified sets' adjusted p-values. The prespecified sets are
# therefore analysed exactly as published, in their own call, and the random
# sets are calibrated separately.
#
# WHAT IT LICENSES
#   LICENSES: a statement that the prespecified enrichments are not a generic
#             consequence of the genome-wide expression shift.
#   DOES NOT: any claim about causality, spatial correspondence with imaging,
#             or that the transcriptomic comparison validates the MRS findings.
#
# INPUTS
#   RES_DIR/deseq2_results_protein_coding.tsv  - DESeq2 output from script 01,
#       one row per tested protein-coding gene, with columns `stat` (Wald
#       statistic) and `symbol` (HGNC). This is the table published as sheet
#       Supp4 of the paper's supplementary workbook.
#   META_DIR/pathways.gmt                      - the 11 prespecified gene sets
#       written by script 00 (five BRETIGEA neuronal tiers, KEGG
#       glycerophospholipid metabolism, Reactome choline metabolism, KEGG
#       glutamatergic synapse, KEGG alanine-aspartate-glutamate metabolism,
#       KEGG oxidative phosphorylation, KEGG TCA cycle).
#
# OUTPUTS
#   RES_DIR/negative_control_gsea.tsv, RES_DIR/negative_control_sessionInfo.txt
#
# Run:  Rscript 04_negative_control_gsea.R
# =============================================================================

suppressMessages({library(fgsea); library(data.table)})
cat("fgsea version:", as.character(packageVersion("fgsea")), "\n\n")

# =============================================================================
# CONFIGURATION
# =============================================================================
USER     <- Sys.getenv("USER")
BASE_DIR <- file.path("/scratch", USER, "glioma_vs_normal")
RES_DIR  <- file.path(BASE_DIR, "results")
META_DIR <- file.path(BASE_DIR, "metadata")

SEED <- 20260814L              # the seed used for the published calibration
set.seed(SEED)
B <- 500L                      # random size-matched sets per prespecified set

RES <- file.path(RES_DIR,  "deseq2_results_protein_coding.tsv")
GMT <- file.path(META_DIR, "pathways.gmt")
OUT <- RES_DIR
for (p in c(RES, GMT)) if (!file.exists(p)) stop("Input not found: ", p)

# ---- ranked list, built exactly as in the main analysis ----------------------
res <- fread(RES)
rk  <- res[!is.na(stat) & !is.na(symbol) & symbol != ""]
rk  <- rk[order(-abs(rk$stat))][!duplicated(symbol)]   # duplicate symbols: keep max |stat|
ranks <- sort(setNames(rk$stat, rk$symbol), decreasing = TRUE)
cat("ranked genes:", length(ranks), "\n")

pw <- gmtPathways(GMT)
cat("prespecified gene sets:", length(pw), "\n\n")

# ---- prespecified sets, analyzed exactly as published ------------------------
obs <- fgsea(pw, ranks, minSize = 5, maxSize = 500, nPermSimple = 10000)
obs <- obs[order(padj)]

# ---- size-matched random negative controls -----------------------------------
universe <- names(ranks)
calib <- rbindlist(lapply(seq_len(nrow(obs)), function(i) {
  nm <- obs$pathway[i]; sz <- obs$size[i]
  rnd <- setNames(lapply(seq_len(B), function(b) sample(universe, sz)),
                  sprintf("RND_%s_%04d", sz, seq_len(B)))
  fr <- suppressWarnings(fgsea(rnd, ranks, minSize = 5, maxSize = 500,
                               nPermSimple = 1000))
  data.table(
    pathway           = nm,
    size              = sz,
    observed_NES      = obs$NES[i],
    observed_padj     = obs$padj[i],
    random_NES_median = median(fr$NES, na.rm = TRUE),
    random_NES_q025   = as.numeric(quantile(fr$NES, .025, na.rm = TRUE)),
    random_NES_q975   = as.numeric(quantile(fr$NES, .975, na.rm = TRUE)),
    random_absNES_max = max(abs(fr$NES), na.rm = TRUE),
    n_random          = nrow(fr),
    calibrated_p      = (sum(abs(fr$NES) >= abs(obs$NES[i]), na.rm = TRUE) + 1) /
                        (nrow(fr) + 1),
    random_frac_FDR05 = mean(fr$padj < 0.05, na.rm = TRUE))
}))

calib[, calibrated_p_BH := p.adjust(calibrated_p, method = "BH")]
calib[, abs_observed_NES := abs(observed_NES)]
calib <- calib[order(calibrated_p, -abs_observed_NES)]

print(calib[, .(pathway, size,
                observed_NES    = round(observed_NES, 2),
                observed_padj   = signif(observed_padj, 3),
                random_NES_95CI = sprintf("%.2f to %.2f", random_NES_q025, random_NES_q975),
                random_absNESmax = round(random_absNES_max, 2),
                calibrated_p    = signif(calibrated_p, 3))])

fwrite(calib, file.path(OUT, "negative_control_gsea.tsv"), sep = "\t")

cat(sprintf("\nB = %d size-matched random sets per prespecified set.\n", B))
cat(sprintf("Prespecified sets whose |NES| exceeded EVERY size-matched random set: %d of %d\n",
            calib[abs(observed_NES) > random_absNES_max, .N], nrow(calib)))
cat(sprintf("Mean proportion of random sets reaching FDR < 0.05 in their own run: %.4f\n",
            mean(calib$random_frac_FDR05, na.rm = TRUE)))

writeLines(capture.output(sessionInfo()), file.path(OUT, "negative_control_sessionInfo.txt"))
cat("\nDone -> ", OUT, "\n")
