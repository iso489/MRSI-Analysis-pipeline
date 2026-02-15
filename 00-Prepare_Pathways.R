#!/usr/bin/env Rscript
# =============================================================================
# Prepare Pathway Gene Sets - OFFLINE VERSION (no internet required)
# Contains hardcoded KEGG/Reactome gene sets + BRETIGEA-derived neuronal markers
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(BRETIGEA)
})

cat("Preparing pathway gene sets for radiogenomic analysis (offline mode)\n\n")

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
BASE_DIR <- "/scratch/ilyaso/tcga_lgg_open"
DATA_DIR <- file.path(BASE_DIR, "data_rnaseq_open")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# BRETIGEA neuronal markers (human-only meta-analysis)
# NOTE:
#  - In BRETIGEA v1.0.3, markers_df_human_brain uses column "cell" (not "cell_type")
#  - We use top-N to preserve specificity (recommended: 200)
#  - Force include NAT8L (NAA synthase)
#  - Exclude ASPA (oligodendrocyte-enriched; NAA catabolism)
# -----------------------------------------------------------------------------
neu_all <- markers_df_human_brain$markers[markers_df_human_brain$cell == "neu"]
neu_all <- unique(na.omit(neu_all))

make_neu <- function(N) {
  neu <- head(neu_all, N)
  neu <- unique(c(neu, "NAT8L"))  # force include
  neu <- setdiff(neu, "ASPA")     # exclude
  neu
}

# Choose the primary neuronal signature here
NEU_PRIMARY_N <- 200
NEURONAL_MARKERS <- make_neu(NEU_PRIMARY_N)

# Optional sensitivity sets (kept in the GMT as additional pathways)
NEURONAL_MARKERS_50  <- make_neu(50)
NEURONAL_MARKERS_100 <- make_neu(100)
NEURONAL_MARKERS_300 <- make_neu(300)
NEURONAL_MARKERS_500 <- make_neu(500)

cat("BRETIGEA neuron markers available (human):", length(neu_all), "\n")
cat("Primary neuronal signature size:", length(NEURONAL_MARKERS), "\n\n")

# -----------------------------------------------------------------------------
# Hardcoded gene sets from KEGG and Reactome
# These are the canonical gene sets for the manuscript hypotheses
# Gene symbols are human (Homo sapiens)
# -----------------------------------------------------------------------------
pathways <- list(

  # =========================================================================
  # 1. MEMBRANE SYNTHESIS AND TURNOVER (Cho-related MRSI signatures)
  # =========================================================================

  # KEGG Glycerophospholipid metabolism (hsa00564)
  KEGG_GLYCEROPHOSPHOLIPID_METABOLISM = c(
    "AGPAT1", "AGPAT2", "AGPAT3", "AGPAT4", "AGPAT5", "AGPAT6", "AGPS",
    "CDIPT", "CEPT1", "CHAT", "CHKA", "CHKB", "CHPT1", "CKI", "CRLS1",
    "DGKA", "DGKB", "DGKD", "DGKE", "DGKG", "DGKH", "DGKI", "DGKK", "DGKQ", "DGKZ",
    "ETNK1", "ETNK2", "GPCPD1", "GPLD1", "GNPAT",
    "LCAT", "LPCAT1", "LPCAT2", "LPCAT3", "LPCAT4", "LPIN1", "LPIN2", "LPIN3",
    "LPGAT1", "MBOAT1", "MBOAT2", "MBOAT7",
    "PCYT1A", "PCYT1B", "PCYT2", "PEMT", "PHOSPHO1", "PISD",
    "PLA2G1B", "PLA2G2A", "PLA2G2C", "PLA2G2D", "PLA2G2E", "PLA2G2F",
    "PLA2G3", "PLA2G4A", "PLA2G4B", "PLA2G4C", "PLA2G4D", "PLA2G4E", "PLA2G4F",
    "PLA2G5", "PLA2G6", "PLA2G10", "PLA2G12A", "PLA2G12B", "PLA2G15",
    "PLCB1", "PLCB2", "PLCB3", "PLCB4", "PLCD1", "PLCD3", "PLCD4",
    "PLCE1", "PLCG1", "PLCG2", "PLCH1", "PLCH2", "PLCZ1",
    "PLD1", "PLD2", "PLD3", "PLD4", "PNPLA6", "PNPLA7", "PNPLA8",
    "PTDSS1", "PTDSS2", "TAZ"
  ),

  # Reactome choline metabolism (your offline set)
  REACTOME_CHOLINE_METABOLISM = c(
    "CHKA", "CHKB", "PCYT1A", "PCYT1B", "CHPT1", "CEPT1",
    "PLD1", "PLD2", "PLA2G4A", "LPCAT1", "LPCAT2",
    "SLC44A1", "SLC44A2", "SLC44A3", "SLC44A4", "SLC44A5",
    "SLC5A7",
    "PEMT", "PCTP", "STARD7", "STARD10",
    "GDPD1", "GDPD2", "GDPD3", "GDPD5", "GPCPD1"
  ),

  # =========================================================================
  # 2. NEURONAL AND AXONAL INTEGRITY (NAA-related MRSI signatures)
  # =========================================================================
  # Primary neuronal identity set derived from BRETIGEA (human-only)
  NEURONAL_MARKERS_BRETIGEA_200 = NEURONAL_MARKERS,

  # Sensitivity pathways
  NEURONAL_MARKERS_BRETIGEA_50  = NEURONAL_MARKERS_50,
  NEURONAL_MARKERS_BRETIGEA_100 = NEURONAL_MARKERS_100,
  NEURONAL_MARKERS_BRETIGEA_300 = NEURONAL_MARKERS_300,
  NEURONAL_MARKERS_BRETIGEA_500 = NEURONAL_MARKERS_500,

  # =========================================================================
  # 3. GLUTAMATE-GLUTAMINE CYCLING (Glu/Cr and Gln/Cr MRSI signatures)
  # =========================================================================

  # KEGG Glutamatergic synapse (hsa04724)
  KEGG_GLUTAMATERGIC_SYNAPSE = c(
    "GRIN1", "GRIN2A", "GRIN2B", "GRIN2C", "GRIN2D", "GRIN3A", "GRIN3B",
    "GRIA1", "GRIA2", "GRIA3", "GRIA4",
    "GRIK1", "GRIK2", "GRIK3", "GRIK4", "GRIK5",
    "GRM1", "GRM2", "GRM3", "GRM4", "GRM5", "GRM6", "GRM7", "GRM8",
    "SLC1A1", "SLC1A2", "SLC1A3", "SLC1A6", "SLC1A7",
    "SLC17A7", "SLC17A6", "SLC17A8",
    "GLUL",
    "GLS", "GLS2",
    "GLUD1", "GLUD2",
    "DLG1", "DLG2", "DLG3", "DLG4",
    "SHANK1", "SHANK2", "SHANK3",
    "HOMER1", "HOMER2", "HOMER3",
    "CAMK2A", "CAMK2B", "CAMK2D", "CAMK2G",
    "PRKCG", "PRKCA", "PRKCB",
    "PLCB1", "PLCB4", "GNAQ", "GNA11"
  ),

  # KEGG Alanine, aspartate and glutamate metabolism (hsa00250)
  KEGG_ALANINE_ASPARTATE_GLUTAMATE_METABOLISM = c(
    "GOT1", "GOT2",
    "GPT", "GPT2",
    "GLS", "GLS2",
    "GLUL",
    "GLUD1", "GLUD2",
    "GAD1", "GAD2",
    "ASNS",
    "ASS1",
    "ASL",
    "ADSL",
    "ADSS1", "ADSS2",
    "ABAT",
    "ALDH5A1",
    "CAD",
    "CPS1",
    "AGXT", "AGXT2",
    "DDO", "DAO"
  ),

  # =========================================================================
  # 4. REDOX AND MITOCHONDRIAL METABOLISM (IDH biological context)
  # =========================================================================

  # KEGG Oxidative phosphorylation (hsa00190)
  KEGG_OXIDATIVE_PHOSPHORYLATION = c(
    "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4", "NDUFA5", "NDUFA6", "NDUFA7",
    "NDUFA8", "NDUFA9", "NDUFA10", "NDUFA11", "NDUFA12", "NDUFA13",
    "NDUFAB1",
    "NDUFB1", "NDUFB2", "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB7",
    "NDUFB8", "NDUFB9", "NDUFB10", "NDUFB11",
    "NDUFC1", "NDUFC2",
    "NDUFS1", "NDUFS2", "NDUFS3", "NDUFS4", "NDUFS5", "NDUFS6", "NDUFS7", "NDUFS8",
    "NDUFV1", "NDUFV2", "NDUFV3",
    "SDHA", "SDHB", "SDHC", "SDHD",
    "UQCRC1", "UQCRC2", "UQCRFS1", "UQCRB", "UQCRH", "UQCRQ",
    "CYC1", "UQCR10", "UQCR11",
    "COX4I1", "COX4I2", "COX5A", "COX5B", "COX6A1", "COX6A2",
    "COX6B1", "COX6B2", "COX6C", "COX7A1", "COX7A2", "COX7B", "COX7C",
    "COX8A", "COX8C",
    "MT-CO1", "MT-CO2", "MT-CO3",
    "ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5F1D", "ATP5F1E",
    "ATP5PB", "ATP5MC1", "ATP5MC2", "ATP5MC3", "ATP5PD", "ATP5ME",
    "ATP5MF", "ATP5MG", "ATP5PF", "ATP5PO", "ATP5IF1",
    "MT-ATP6", "MT-ATP8",
    "CYCS"
  ),

  # KEGG Citrate cycle / TCA cycle (hsa00020)
  KEGG_CITRATE_CYCLE_TCA_CYCLE = c(
    "PC",
    "PDHA1", "PDHA2", "PDHB",
    "DLAT", "DLD",
    "PDK1", "PDK2", "PDK3", "PDK4",
    "PDP1", "PDP2",
    "CS",
    "ACO1", "ACO2",
    "IDH1", "IDH2", "IDH3A", "IDH3B", "IDH3G",
    "OGDH", "OGDHL",
    "DLST",
    "SUCLG1", "SUCLG2", "SUCLA2",
    "SDHA", "SDHB", "SDHC", "SDHD",
    "FH",
    "MDH1", "MDH2",
    "PCK1", "PCK2",
    "ME1", "ME2", "ME3",
    "ACLY"
  )
)

# -----------------------------------------------------------------------------
# Validate and deduplicate gene sets
# -----------------------------------------------------------------------------
cat("Validating gene sets:\n")
for (nm in names(pathways)) {
  pathways[[nm]] <- unique(pathways[[nm]])
  cat(sprintf("  %s: %d genes\n", nm, length(pathways[[nm]])))
}

# -----------------------------------------------------------------------------
# Write GMT file
# GMT format: pathway_name <TAB> description <TAB> gene1 <TAB> gene2 <TAB> ...
# -----------------------------------------------------------------------------
gmt_file <- file.path(DATA_DIR, "pathways.gmt")
cat("\nWriting GMT file:", gmt_file, "\n")

con <- file(gmt_file, "w")
for (pathway_name in names(pathways)) {
  genes <- pathways[[pathway_name]]
  line <- paste(c(pathway_name, pathway_name, genes), collapse = "\t")
  writeLines(line, con)
}
close(con)

cat("GMT file created successfully!\n")

# Verify
cat("\nVerifying GMT file:\n")
gmt_lines <- readLines(gmt_file)
for (line in gmt_lines) {
  parts <- strsplit(line, "\t")[[1]]
  cat(sprintf("  %s: %d genes\n", parts[1], length(parts) - 2))
}

cat("\n=== Pathway preparation complete ===\n")
