#!/bin/bash
#SBATCH --job-name=diagnoseHinv
#SBATCH --output=diagnoseHinv.out
#SBATCH --error=diagnoseHinv.err
#SBATCH --partition=GPU            # same high-memory partition
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=100G
#SBATCH --time=12:00:00

set -euo pipefail

# -----------------------------
# Load Conda environment
# -----------------------------
source /install/software/2024/minpy3/python3.10/etc/profile.d/conda.sh
conda activate /data/Genetics/analysis/R1681_OviSeq/Dermot/envs/dermo_python_env

# -----------------------------
# Paths
# -----------------------------
DIR="$HOME/Dermot_analysis/Phd/Paper_3"
SCRIPT="$DIR/scripts/diagnostic_hinv.py"
HINV_FILE="$DIR/scripts/sheep50k_fast_Hinv_PD.giv"

mkdir -p "$DIR/results"

# -----------------------------
# Parallelization environment
# -----------------------------
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export BLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export PYTHONUNBUFFERED=1

# -----------------------------
# Run
# -----------------------------
echo "=== Starting diagnostic_hinv.py ==="
date
echo "Hinv file:      $HINV_FILE"
echo "Threads:        $SLURM_CPUS_PER_TASK"
echo "Memory:         $SLURM_MEM_PER_NODE MB"
echo "==================================="

srun python3 -u "$SCRIPT" "$HINV_FILE"

echo "=== H^-1 diagnostic finished $(date) ==="
