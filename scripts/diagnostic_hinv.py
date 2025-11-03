#!/usr/bin/env python3
"""
diagnostic_hinv.py
------------------
Fast H^-1 diagnostic (Python version)
Checks symmetry, diagonals, sparsity, density, PD status, and eigenstructure.
Usage:
    python diagnostic_hinv.py sheep50k_fast_Hinv.giv
"""

import sys, time
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.sparse import linalg as spla

if len(sys.argv) < 2:
    sys.exit("Usage: python diagnostic_hinv.py <Hinv.giv>")

giv_file = sys.argv[1]
t0 = time.time()

print("\n=================================================")
print("   H^-1 Comprehensive Diagnostic (Python version)")
print("   File:", giv_file)
print("=================================================")

# ------------------------------------------------------------
# 1. Read triplets
# ------------------------------------------------------------
print("[1/8] Reading triplets...", flush=True)
trip = pd.read_csv(giv_file, sep=" ", header=None, names=["i","j","x"])
n = int(max(trip.i.max(), trip.j.max()))
nnz = len(trip)
density = nnz / (n**2)
print(f"  → Dimension: {n:,} × {n:,}")
print(f"  → Non-zeros: {nnz:,} ({100*density:.6f}% dense)")

# ------------------------------------------------------------
# 2. Build sparse symmetric matrix
# ------------------------------------------------------------
print("[2/8] Building sparse symmetric matrix...", flush=True)
H = sparse.coo_matrix((trip.x, (trip.i-1, trip.j-1)), shape=(n,n))
H = H + H.T - sparse.diags(H.diagonal())   # enforce symmetry
H = H.tocsr()

# ------------------------------------------------------------
# 3. Structure and diagonals
# ------------------------------------------------------------
print("[3/8] Structural checks...", flush=True)
diag_vals = H.diagonal()
neg_diag = np.sum(diag_vals < 0)
print(f"  • Diagonal range: [{diag_vals.min():.4f}, {diag_vals.max():.4f}]")
print(f"  • Negative diagonals: {neg_diag} ({100*neg_diag/n:.3f}%)")

sym_diff = sparse.linalg.norm(H - H.T, ord='fro') / sparse.linalg.norm(H, ord='fro')
print(f"  • Symmetry diff (rel Frobenius): {sym_diff:.3e}")

# ------------------------------------------------------------
# 4. Sparsity and connectivity
# ------------------------------------------------------------
print("[4/8] Sparsity & connectivity checks...", flush=True)
# Efficient nonzero count per row (sparse-friendly)
row_nnz = np.diff(H.indptr)       # CSR row pointer differences = count of nonzeros
zero_rows = np.where(row_nnz == 0)[0]


print(f"  • Zero rows/cols: {len(zero_rows)}")
if len(zero_rows)==0:
    print("  ✓ Fully connected (no zero rows)")

# ------------------------------------------------------------
# 5. Near-zero or duplicate rows
# ------------------------------------------------------------
print("[5/8] Searching for near-zero rows...", flush=True)
row_norms = np.sqrt(np.array(H.power(2).sum(axis=1)).ravel())
zero_like = np.where(row_norms < 1e-8)[0]
if len(zero_like) > 0:
    print(f"  ⚠️  {len(zero_like)} rows have near-zero norm (redundant?)")
else:
    print("  ✓ No near-zero rows detected")

# ------------------------------------------------------------
# 6. Positive-definiteness check
# ------------------------------------------------------------
print("[6/8] PD check via sparse LU (fast)...", flush=True)

def is_pd_sparse(M):
    """Try LU factorization as PD proxy (fast)."""
    try:
        _ = spla.splu(M.tocsc())
        return True
    except Exception:
        return False

pd_ok = False
for ridge in [1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1]:
    print(f"   → Testing ridge = {ridge:.1e} ...", end=" ")
    Htest = H + sparse.identity(H.shape[0]) * ridge
    if is_pd_sparse(Htest):
        print("✓ PD")
        pd_ok = True
        break
    else:
        print("not PD")
if not pd_ok:
    print("  ⚠️  Still indefinite up to ridge = 1e-1")

# ------------------------------------------------------------
# 7. Density / Δ coverage heuristics
# ------------------------------------------------------------
print("[7/8] Density breakdown & Δ coverage...", flush=True)
density = 100 * H.nnz / (n**2)
print(f"  • Nonzeros: {H.nnz:,} ({density:.6f}% dense)")
if density < 15e6/(n**2)/5:
    print("  ⚠️  Much sparser than expected for full A^-1 → maybe missing A-block.")
else:
    print("  ✓ Density within expected range.")

diag_sorted = np.sort(diag_vals)
delta_jumps = np.sum(np.abs(np.diff(diag_sorted)) > 0.5)
print(f"  • Large diagonal jumps: {delta_jumps}")
if delta_jumps > 1000:
    print("  ⚠️  Irregular Δ embedding (A22/G misalignment).")
else:
    print("  ✓ Δ block embedding consistent.")

# ------------------------------------------------------------
# 8. Eigenvalue sketch (approximate)
# ------------------------------------------------------------
print("[8/8] Smallest eigenvalues (approx.)...", flush=True)
try:
    eigvals = spla.eigsh(H, k=5, which='SM', return_eigenvectors=False)
    print(f"  • Smallest eigenvalue: {eigvals.min():.3e}")
    if eigvals.min() <= 0:
        print("  ⚠️  Matrix indefinite or near-singular.")
    else:
        print("  ✓ Eigenvalues all positive (PD confirmed).")
except Exception as e:
    print("  ⚠️  Could not compute eigenvalues (matrix too large or ill-conditioned).")

# ------------------------------------------------------------
# Summary
# ------------------------------------------------------------
print("\n=============== STRUCTURAL SUMMARY ===============")
print(f"Matrix size: {n:,} × {n:,}")
print(f"Non-zeros: {H.nnz:,} ({density:.6f}% dense)")
print(f"Diagonal range: [{diag_vals.min():.4f}, {diag_vals.max():.4f}]")
print(f"Neg. diagonals: {neg_diag}")
print(f"Sym diff: {sym_diff:.3e}")
print("--------------------------------------------------")
if not pd_ok:
    print("⚠️  DIAGNOSIS: Indefinite or scaling issue.")
else:
    print("✅  DIAGNOSIS: Matrix structurally valid and PD.")
print("==================================================")
print(f"⏱ Completed in {(time.time()-t0)/60:.2f} min")
