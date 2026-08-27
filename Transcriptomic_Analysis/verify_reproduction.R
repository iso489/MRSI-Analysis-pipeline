#!/usr/bin/env Rscript
# =============================================================================
# verify_reproduction.R - automated check of a fresh run against EXPECTED_RESULTS
# =============================================================================
#
# Run this after 00-05. It reads the outputs a run just produced, compares them
# with EXPECTED_RESULTS.tsv, prints a PASS/FAIL line per quantity and exits
# non-zero on any failure, so it can gate CI.
#
# Enrichment scores come from a permutation procedure and are reproducible only
# to within Monte-Carlo error, so each quantity carries its own tolerance;
# integer counts carry a tolerance of 0 unless the DESeq2 fit is expected to
# wobble at the margin of significance.
#
# Run:  Rscript verify_reproduction.R
# =============================================================================

suppressMessages(library(data.table))

HERE <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
USER <- Sys.getenv("USER")
BASE <- Sys.getenv("GVN_BASE", unset = file.path("/scratch", USER, "glioma_vs_normal"))
RES  <- file.path(BASE, "results")
DON  <- file.path(RES, "donor_level")

exp_file <- file.path(HERE, "EXPECTED_RESULTS.tsv")
if (!file.exists(exp_file)) stop("missing EXPECTED_RESULTS.tsv beside this script")
expected <- fread(exp_file)

get <- function(q) {
  if (q == "gtex_samples")
    return(fread(file.path(HERE, "sample_donor_manifest.tsv"))[condition == "normal", .N])
  if (q == "gtex_donors")
    return(uniqueN(fread(file.path(HERE, "sample_donor_manifest.tsv"))[condition == "normal", donor]))
  if (q == "gtex_donors_with_repeats") {
    m <- fread(file.path(HERE, "sample_donor_manifest.tsv"))[condition == "normal"]
    return(m[, .N, by = donor][N > 1, .N])
  }
  if (q == "tcga_participants")
    return(uniqueN(fread(file.path(HERE, "sample_donor_manifest.tsv"))[condition == "glioma", donor]))
  if (q == "protein_coding_genes_tested")
    return(nrow(fread(file.path(RES, "deseq2_results_protein_coding.tsv"))))
  if (startsWith(q, "n_DE_")) {
    d <- fread(file.path(DON, "de_summary.tsv"))
    key <- c(sample_level = "A0 sample-level", donor_level = "A1 donor-level (PRIMARY)",
             one_per_donor = "A2 one-per-donor")[sub("^n_DE_", "", q)]
    return(d[analysis == key, n_DE])
  }
  if (startsWith(q, "NES_")) {
    f <- fread(file.path(DON, "fgsea_all_analyses.tsv"))
    return(f[analysis == "A1 donor-level (PRIMARY)" & pathway == sub("^NES_", "", q), NES])
  }
  if (q == "NAT8L_log2FC_donor_level") {
    s <- fread(file.path(DON, "single_genes_all_analyses.tsv"))
    return(s[analysis == "A1 donor-level (PRIMARY)" & symbol == "NAT8L", log2FoldChange])
  }
  if (q == "sets_exceeding_all_random_sets")
    return(fread(file.path(DON, "A1_negative_control_gsea.tsv"))[calibrated_p <= 0.002, .N])
  NA_real_
}

cat(sprintf("%-46s %12s %12s %8s  %s\n", "QUANTITY", "EXPECTED", "OBSERVED", "TOL", "RESULT"))
fails <- 0L
for (i in seq_len(nrow(expected))) {
  q   <- expected$quantity[i]
  exp <- as.numeric(expected$expected[i])
  tol <- as.numeric(expected$tolerance[i])
  obs <- tryCatch(as.numeric(get(q)), error = function(e) NA_real_)
  ok  <- !is.na(obs) && abs(obs - exp) <= tol
  if (!ok) fails <- fails + 1L
  cat(sprintf("%-46s %12s %12s %8s  %s\n", q, format(exp), 
              if (is.na(obs)) "MISSING" else format(round(obs, 3)),
              format(tol), if (ok) "PASS" else "FAIL"))
}
cat(sprintf("\n%d quantities checked, %d failures\n", nrow(expected), fails))
if (fails > 0L) quit(status = 1L)
cat("Reproduction verified against EXPECTED_RESULTS.tsv\n")
