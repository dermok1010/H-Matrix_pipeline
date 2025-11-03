#!/bin/bash
#SBATCH --job-name=makeGinv
#SBATCH --output=makeGinv_%j.out
#SBATCH --error=makeGinv_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dermot.kelly@teagasc.ie

set -euo pipefail
module load apptainer/1.1.9

# ============================================================
# CONFIGURATION
# ============================================================
BASE="$HOME/Dermot_analysis/Phd/Paper_3/H_matrix"
SCRIPT="$BASE/scripts/build_ginv.R"
IMG="$HOME/Dermot_analysis/tidyverse_4.2.sif"

# Input data
PED="$BASE/data/pedigree_full_sas_style.csv"
GRM_PREFIX="$BASE/data/sheep50k"             # expects sheep50k.grm + sheep50k.grm.id
OUT_PREFIX="$BASE/results/sheep50k_Ginv"     # output .giv + .order

# Bind directories for Singularity container
BIND_OPTS="--bind $HOME --bind /data/Genetics"

# Output directories
mkdir -p "$BASE/results" "$BASE/logs"

# Thread control
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export BLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

# ============================================================
# RUN
# ============================================================
echo "=== Starting G^-1 build job on $(date) ==="
echo "Script:    $SCRIPT"
echo "Pedigree:  $PED"
echo "GRM:       $GRM_PREFIX"
echo "Output:    $OUT_PREFIX"
echo "Threads:   $SLURM_CPUS_PER_TASK"
echo "Memory:    $SLURM_MEM_PER_NODE MB"
echo "==========================================================="

srun apptainer exec \
  --env OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK \
  --env MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK \
  --env OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK \
  --env BLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK \
  $BIND_OPTS "$IMG" \
  Rscript "$SCRIPT" "$PED" "$GRM_PREFIX" "$OUT_PREFIX"

echo "=== ✅ G^-1 build finished successfully on $(date) ==="
