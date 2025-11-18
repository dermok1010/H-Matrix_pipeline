#!/bin/bash
#SBATCH --job-name=fixGinv
#SBATCH --output=fixGinv.out
#SBATCH --error=fixGinv.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
#SBATCH --time=12:00:00

set -euo pipefail
module load apptainer/1.1.9

# ============================================================
# CONFIGURATION
# ============================================================
BASE="$HOME/Dermot_analysis/Phd/"
SCRIPT="$BASE/Noirin_H/fix_g.R"
IMG="$HOME/Dermot_analysis/tidyverse_4.2.sif"

# Bind directories for Apptainer container
BIND_OPTS="--bind $HOME --bind /data/Genetics"

# Thread control
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export BLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

# ============================================================
# RUN
# ============================================================
echo "=== Starting G^-1 FIX job on $(date) ==="
echo "Script:     $SCRIPT"
echo "Container:  $IMG"
echo "Threads:    $SLURM_CPUS_PER_TASK"
echo "Memory:     $SLURM_MEM_PER_NODE MB"
echo "==========================================================="

srun apptainer exec \
  --env OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK \
  --env MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK \
  --env OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK \
  --env BLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK \
  $BIND_OPTS "$IMG" \
  Rscript "$SCRIPT"

echo "==========================================================="
echo "✅ Job completed successfully on $(date)"
echo "==========================================================="
