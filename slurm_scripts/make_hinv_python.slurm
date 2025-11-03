#!/bin/bash
#SBATCH --job-name=makeHinvPy
#SBATCH --output=makeHinvPy.out
#SBATCH --error=makeHinvPy.err
#SBATCH --partition=GPU           # high-memory node
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=300G
#SBATCH --time=24:00:00

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
SCRIPT="$DIR/genotypes/lets_make_Hinv_py.py"

AINV="$DIR/results/sheep50k_Ainv.mtx"
A22INV="$DIR/results/sheep50k_A22inv.mtx"
GINV="$DIR/scripts/Ginv_sheep50k_SPARSE.giv"
GORDER="$DIR/scripts/Ginv_sheep50k_SPARSE.order"
OUTPREFIX="$DIR/results/sheep50k_fast"

mkdir -p "$DIR/results"

# -----------------------------
# Parallelization environment
# -----------------------------
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export BLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export PYTHONUNBUFFERED=1  # ensure Python prints live

# -----------------------------
# Run
# -----------------------------
echo "=== Starting lets_make_Hinv_py.py ==="
date
echo "Ainv file:      $AINV"
echo "A22inv file:    $A22INV"
echo "Ginv file:      $GINV"
echo "Order file:     $GORDER"
echo "Output prefix:  $OUTPREFIX"
echo "Threads:        $SLURM_CPUS_PER_TASK"
echo "Memory:         $SLURM_MEM_PER_NODE MB"
echo "==================================="

srun python3 -u "$SCRIPT" \
  "$AINV" \
  "$A22INV" \
  "$GINV" \
  "$GORDER" \
  "$OUTPREFIX"

echo "=== H^-1 build finished $(date) ==="
