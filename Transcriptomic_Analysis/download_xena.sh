#!/bin/bash
# =============================================================================
# 00_download_xena.sh
# Downloads Toil/CGL reprocessed expression data from UCSC Xena
# Run this on a login node (or transfer node) BEFORE submitting Analysis.sbatch
# =============================================================================

set -euo pipefail

DL_DIR="/scratch/${USER}/glioma_vs_normal"
mkdir -p "$DL_DIR"

XENA_BASE="https://toil-xena-hub.s3.us-east-1.amazonaws.com/download"

FILES=(
  "TcgaTargetGtex_gene_expected_count.gz"
  "TcgaTargetGTEX_phenotype.txt.gz"
)

for f in "${FILES[@]}"; do
  dest="${DL_DIR}/${f}"

  # If the file exists, verify integrity first
  if [[ -f "$dest" ]]; then
    echo "[CHECK] Verifying $f ..."
    if gzip -t "$dest" 2>/dev/null; then
      echo "[OK]    $f is intact — skipping download."
      continue
    else
      echo "[CORRUPT] $f failed gzip integrity check — re-downloading."
      rm -f "$dest"
    fi
  fi

  echo "[DOWNLOADING] $f ..."
  curl -L --retry 5 --retry-delay 10 -C - -o "$dest" "${XENA_BASE}/${f}"

  # Verify after download
  echo "[CHECK] Verifying $f after download ..."
  if gzip -t "$dest" 2>/dev/null; then
    echo "[OK]    $f downloaded and verified."
  else
    echo ""
    echo "[ERROR] $f is CORRUPT after download."
    echo "        This usually means the connection was interrupted."
    echo "        Delete the file and try again, or use a transfer node:"
    echo "          rm $dest"
    echo "          ssh ${USER}@narval-dtn  # or beluga-dtn, cedar-dtn, etc."
    echo "          bash 00_download_xena.sh"
    rm -f "$dest"
    exit 1
  fi
done

echo ""
echo "All files saved to: $DL_DIR"
ls -lh "$DL_DIR"/*.gz
echo ""
echo "You can now submit: sbatch Analysis.sbatch"
