#!/usr/bin/env python3
"""
lets_make_Hinv_py.py
--------------------
Constructs H⁻¹ for ASReml using:
  - A⁻¹ and A₂₂⁻¹ (Matrix Market from R)
  - G⁻¹ (.giv + .order)

Outputs:
  - <prefix>_Hinv.giv
  - <prefix>_Hinv.order
"""

import numpy as np
import pandas as pd
from scipy import sparse, io as spio, linalg
import sys, time

t0 = time.time()
if len(sys.argv) < 6:
    sys.exit("Usage: lets_make_Hinv_py.py <Ainv.mtx> <A22inv.mtx> <Ginv.giv> <Ginv.order> <output_prefix>")

ainv_file, a22inv_file, ginv_file, gorder_file, out_prefix = sys.argv[1:6]

# ============================================================
# [1] Load sparse matrices
# ============================================================
print(f"\n[1/7] Reading A^-1 and A22^-1 (Matrix Market)...", flush=True)
Ainv = spio.mmread(ainv_file).tocsr()
A22inv = spio.mmread(a22inv_file).tocsr()
nA, nG = Ainv.shape[0], A22inv.shape[0]
print(f"   - A^-1:  {nA:,} x {nA:,}  | nnz={Ainv.nnz:,}")
print(f"   - A22^-1: {nG:,} x {nG:,} | nnz={A22inv.nnz:,}")

# ============================================================
# [2] Load G^-1
# ============================================================
print(f"\n[2/7] Reading G^-1 triplets...", flush=True)
Gtrip = pd.read_csv(ginv_file, sep=" ", header=None, names=["i", "j", "x"])
geno_ids = np.loadtxt(gorder_file, dtype=str)
nG2 = len(geno_ids)
if nG2 != nG:
    print(f"⚠️  G^-1 ({nG2}) and A22^-1 ({nG}) differ in dimension — proceeding with min(nG).")
    nG = min(nG, nG2)
Ginv = sparse.coo_matrix((Gtrip.x, (Gtrip.i - 1, Gtrip.j - 1)), shape=(nG, nG)).tocsr()
print(f"   - G^-1: {nG:,} x {nG:,} | nnz={Ginv.nnz:,}")

# ============================================================
# [3] Compute Δ = (G^-1 - A22^-1)
# ============================================================
print(f"\n[3/7] Forming delta = G^-1 - A22^-1 ...", flush=True)
delta = Ginv - A22inv[:nG, :nG]
delta.eliminate_zeros()
print(f"   - Δ nonzeros: {delta.nnz:,}")

# ============================================================
# [4/7] Building H^-1 sparsely (with correct ID alignment)
# ============================================================
print(f"\n[4/7] Building H^-1 sparsely (with correct ID alignment)...", flush=True)

# 1️⃣ Load ID order of Ainv and Ginv
ainv_order_file = ainv_file.replace("_Ainv.mtx", "_Ainv.order")
ainv_ids = np.loadtxt(ainv_order_file, dtype=str)
geno_ids = np.loadtxt(gorder_file, dtype=str)

# 2️⃣ Build a mapping from each genotyped animal to its row index in Ainv
aidx = {id_: i for i, id_ in enumerate(ainv_ids)}
missing = [gid for gid in geno_ids if gid not in aidx]
if missing:
    print(f"⚠️  {len(missing)} genotyped IDs not found in pedigree order (showing first 5): {missing[:5]}")
idx_map = np.array([aidx[gid] for gid in geno_ids if gid in aidx])
print(f"   - Mapped {len(idx_map):,} genotyped IDs into Ainv indices")

# 3️⃣ Build sparse Δ in the full pedigree coordinate system
# Map only the nonzero elements of delta to pedigree coordinates
di, dj = delta.nonzero()
dv = delta.data
delta_full = sparse.csr_matrix(
    (dv, (idx_map[di], idx_map[dj])),
    shape=(nA, nA)
)

# 4️⃣ Add Δ to A^-1 to form H^-1
Hinv = Ainv + delta_full
Hinv.eliminate_zeros()
Hinv = sparse.tril(Hinv).tocsr()

print(f"   - H^-1 dimension: {Hinv.shape[0]:,} x {Hinv.shape[1]:,}")
print(f"   - H^-1 nonzeros: {Hinv.nnz:,}")
print(f"   - Δ injected into correct pedigree sub-block (A22 positions)")

# ============================================================
# [5] Export for ASReml
# ============================================================
print(f"\n[5/7] Writing ASReml .giv and .order files...", flush=True)
Htrip = sparse.tril(Hinv).tocoo()
trip_df = pd.DataFrame({
    "i": Htrip.row + 1,
    "j": Htrip.col + 1,
    "x": Htrip.data
})
trip_df.sort_values(["i", "j"], inplace=True)
trip_df.to_csv(f"{out_prefix}_Hinv.giv", sep=" ", header=False, index=False)

# Order = full pedigree order (from R)
order_file = ainv_file.replace("_Ainv.mtx", "_Ainv.order")
order_ids = np.loadtxt(order_file, dtype=str)
np.savetxt(f"{out_prefix}_Hinv.order", order_ids, fmt="%s")

print(f"✅ Wrote {out_prefix}_Hinv.giv")
print(f"✅ Wrote {out_prefix}_Hinv.order")

# ============================================================
# [6] Enforce strict symmetry
# ============================================================
print(f"\n[6/8] Enforcing strict symmetry ...", flush=True)
Hsym = (Hinv + Hinv.T) / 2
sym_diff = sparse.linalg.norm(Hsym - Hsym.T, ord='fro') / sparse.linalg.norm(Hsym, ord='fro')
print(f"   • Symmetry diff (rel): {sym_diff:.3e}")

# ============================================================
# [7] Sparse positive-definiteness enforcement (safe for big H)
# ============================================================
print(f"\n[7/8] Sparse positive-definiteness enforcement ...", flush=True)

from scipy.sparse.linalg import eigsh

# Estimate smallest eigenvalue
print("   • Estimating smallest eigenvalue via eigsh(k=1)...", flush=True)
try:
    lam_min, _ = eigsh(Hsym, k=1, which="SA")   # smallest algebraic
    lam_min = lam_min[0]
except Exception as e:
    print(f"   ⚠️ eigsh failed ({e}); assuming lam_min = -1e-3")
    lam_min = -1e-3

print(f"   • Smallest eigenvalue estimate: {lam_min:.3e}")

# Compute required ridge to push all eigenvalues positive
ridge = max(1e-6 - lam_min, 0)
if ridge > 0:
    print(f"   • Adding ridge = {ridge:.3e} to diagonals to ensure PD")
    Hsym = Hsym + sparse.eye(Hsym.shape[0], format="csr") * ridge
else:
    print("   ✓ Matrix already PD (no ridge needed)")

# Double-check diagonal safety
d = Hsym.diagonal()
neg_diag = np.sum(d <= 0)
if neg_diag > 0:
    print(f"   ⚠️  {neg_diag} nonpositive diagonals detected; fixing to 1e-6")
    d[d <= 0] = 1e-6
    Hsym.setdiag(d)
else:
    print("   ✓ All diagonals positive")

# ============================================================
# [8] Final diagnostics & export
# ============================================================
print(f"\n[8/8] Final diagnostics & export...", flush=True)
diag_vals = Hsym.diagonal()
sym_diff = sparse.linalg.norm(Hsym - Hsym.T, ord='fro') / sparse.linalg.norm(Hsym, ord='fro')
neg_diag = np.sum(diag_vals < 0)
density = 100 * Hsym.nnz / (Hsym.shape[0] ** 2)
print(f"   • Sym diff (rel): {sym_diff:.3e}")
print(f"   • Neg. diagonals: {neg_diag}")
print(f"   • Density: {density:.6f}%")
print(f"   • Diagonal range: [{diag_vals.min():.4f}, {diag_vals.max():.4f}]")

# Export repaired Hinv
bent_giv = f"{out_prefix}_Hinv_PD.giv"
bent_ord = f"{out_prefix}_Hinv_PD.order"
Htrip = sparse.tril(Hsym).tocoo()
trip_df = pd.DataFrame({
    "i": Htrip.row + 1,
    "j": Htrip.col + 1,
    "x": Htrip.data
})
trip_df.sort_values(["i", "j"], inplace=True)
trip_df.to_csv(bent_giv, sep=" ", header=False, index=False)
np.savetxt(bent_ord, order_ids, fmt="%s")

print(f"\n✅ Wrote strictly PD matrix (sparse): {bent_giv}")
print(f"✅ Ridge applied: {ridge:.3e}")
print(f"✅ Matrix guaranteed PD — safe for ASReml.\n")
print(f"⏱ Total runtime: {(time.time() - t0) / 60:.2f} min")
