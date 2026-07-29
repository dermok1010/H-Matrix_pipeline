#!/bin/bash
#SBATCH --job-name=makeGinv
#SBATCH --output=makeGinv_upd.out
#SBATCH --error=makeGinv_upd.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=120G
#SBATCH --time=24:00:00

set -euo pipefail
module load apptainer/1.1.9

# ============================================================
# CONFIGURATION
# ============================================================
BASE="$HOME/Dermot_analysis/Phd/Paper_3/H_matrix"
SCRIPT="$BASE/scripts/build_ginv.R"
IMG="$HOME/Dermot_analysis/tidyverse_4.2.sif"

# Input data
PED="/home/dermot.kelly/Dermot_analysis/Phd/Paper_3/trial/upd_ped_fixed.csv.SRT"
GRM_PREFIX="/home/dermot.kelly/Dermot_primary/refined_paper_3/data/sheep50k_p2"             # expects sheep50k.grm + sheep50k.grm.id
OUT_PREFIX="/home/dermot.kelly/Dermot_analysis/Phd/Paper_3/trial/upd_sheep50k_Ginv"

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
