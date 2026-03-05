#!/bin/bash
#SBATCH --job-name=RNA_seq
#SBATCH --account=[enter_account_name]
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=32G
#SBATCH --output=[enter_output_path]_%j.out
#SBATCH --error=[enter_output_path]_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=[enter_email]

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
module load StdEnv/2023 r/4.5.0

# Prevent OpenBLAS/OpenMP from oversubscribing — let R handle parallelism
export OMP_NUM_THREADS=1

mkdir -p /scratch/ilyaso/glioma_vs_normal/logs

SCRIPTS=/scratch/ilyaso/glioma_vs_normal/scripts

echo "========================================================"
echo " Job:    $SLURM_JOB_ID"
echo " Node:   $SLURMD_NODENAME"
echo " CPUs:   $SLURM_CPUS_PER_TASK"
echo " Mem:    $SLURM_MEM_PER_NODE MB"
echo " Start:  $(date)"
echo "========================================================"

# ---------------------------------------------------------------------------
# 00 Prepare pathway gene sets
# ---------------------------------------------------------------------------
echo ""
echo "--- Script 00: Prepare pathways ---"
Rscript "${SCRIPTS}/00_prepare_pathways.R"
if [ $? -ne 0 ]; then echo "ERROR in script 00"; exit 1; fi

# ---------------------------------------------------------------------------
# 01 DESeq2 glioma vs normal
# ---------------------------------------------------------------------------
echo ""
echo "--- Script 01: DESeq2 ---"
Rscript "${SCRIPTS}/01_deseq2_glioma_vs_normal.R"
if [ $? -ne 0 ]; then echo "ERROR in script 01"; exit 1; fi

# ---------------------------------------------------------------------------
# 02 fgsea enrichment 
# ---------------------------------------------------------------------------
echo ""
echo "--- Script 02: fgsea ---"
Rscript "${SCRIPTS}/02_fgsea_analysis.R"
if [ $? -ne 0 ]; then echo "ERROR in script 02"; exit 1; fi

# ---------------------------------------------------------------------------
# 03 GSVA per-sample pathway activity 
# ---------------------------------------------------------------------------
echo ""
echo "--- Script 03: GSVA ---"
Rscript "${SCRIPTS}/03_gsva_analysis.R"
if [ $? -ne 0 ]; then echo "ERROR in script 03"; exit 1; fi

echo ""
echo "========================================================"
echo " End:    $(date)"
echo " All scripts completed successfully"
echo "========================================================"
