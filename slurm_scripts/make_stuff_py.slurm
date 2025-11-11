#!/bin/bash
#SBATCH --job-name=exportAinv
#SBATCH --output=exportAinv.out
#SBATCH --error=exportAinv.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16       # 32 threads is plenty for makeAinv()
#SBATCH --mem=64G               # ~200 GB safe for large pedigree
#SBATCH --time=12:00:00

set -euo pipefail
module load apptainer/1.1.9

# -------------------------------
# Paths
# -------------------------------
DIR="$HOME/Dermot_analysis/Phd/"
SCRIPT="$DIR/Paper_3/H_matrix/scripts/make_stuff_py.R"     # new lightweight R script
IMG="$HOME/Dermot_analysis/tidyverse_4.2.sif"

PED_FILE="$DIR/Noirin_H/pedigree_corrected.csv"
GENO_IDS="$DIR/Noirin_H/sheep_pedigree_geno.grm.id"
OUTDIR="$DIR/py_results"
OUTPREFIX="$OUTDIR/sheep50k"                 # base prefix for output files

mkdir -p "$OUTDIR"

# -------------------------------
# Threads and BLAS control
# -------------------------------
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK
export BLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK

BIND_OPTS="--bind $HOME --bind /data/Genetics"

# -------------------------------
# Run job
# -------------------------------
echo "=== Starting export_for_python.R job $(date) ==="
echo "Pedigree:      $PED_FILE"
echo "Genotyped IDs: $GENO_IDS"
echo "Output prefix: $OUTPREFIX"
echo "Threads:       $SLURM_CPUS_PER_TASK"
echo "Memory:        $SLURM_MEM_PER_NODE MB"
echo "==============================================="

srun apptainer exec \
  --env OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK \
  --env MKL_NUM_THREADS=$SLURM_CPUS_PER_TASK \
  --env OPENBLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK \
  --env BLAS_NUM_THREADS=$SLURM_CPUS_PER_TASK \
  $BIND_OPTS "$IMG" Rscript "$SCRIPT" \
  "$PED_FILE" "$GENO_IDS" "$OUTPREFIX"

echo "=== Export complete $(date) ==="
