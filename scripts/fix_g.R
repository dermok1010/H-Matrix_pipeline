#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(nadiv)
})
setwd("/home/dermot.kelly/Dermot_analysis/Phd/Paper_3/scripts")
# ============================================================
# CONFIG (adjust only if needed)
# ============================================================
giv_file    <- "Ginv_sheep50k.giv"         # Existing G^-1
order_file  <- "Ginv_sheep50k.order"       # Existing order
ped_file    <- "pedigree_full_sas_style.csv"  # Pedigree

out_prefix  <- "Ginv_sheep50k_FIXED"       # NEW output prefix
out_giv     <- paste0(out_prefix, ".giv")
out_order   <- paste0(out_prefix, ".order")

cat("\n=====================================================\n")
cat("  FIXING EXISTING G^-1 FOR ASREML\n")
cat("=====================================================\n")

cat("Input Ginv:      ", giv_file, "\n")
cat("Order file:      ", order_file, "\n")
cat("Pedigree file:   ", ped_file, "\n")
cat("Output prefix:   ", out_prefix, "\n")
cat("=====================================================\n")

# ============================================================
# 1) Load order (IDs)
# ============================================================
cat("\n[1/6] Reading order file...\n")
ids <- readLines(order_file)
n <- length(ids)
cat("   - Genotyped animals:", n, "\n")

# ============================================================
# 2) Load sparse Ginv from .giv
# ============================================================
cat("\n[2/6] Reading Ginv triplets...\n")
Gtrip <- fread(giv_file, header=FALSE)
if (ncol(Gtrip) < 3) stop("Ginv file must have 3 columns: i j value")
setnames(Gtrip, c("i","j","value"))
cat("   - Non-zero entries in original Ginv:", nrow(Gtrip), "\n")

cat("   - Converting to sparse matrix...\n")
Ginv_sp <- sparseMatrix(
  i = Gtrip$i,
  j = Gtrip$j,
  x = Gtrip$value,
  dims = c(n,n)
)
rownames(Ginv_sp) <- ids
colnames(Ginv_sp) <- ids
rm(Gtrip); gc()

cat("   - Converting to dense matrix (~13GB)...\n")
Ginv <- as.matrix(Ginv_sp)
rm(Ginv_sp); gc()

# Check basic stats
diagG <- diag(Ginv)
cat("   - mean diag(Ginv) BEFORE fix:", mean(diagG), "\n")

# ============================================================
# 3) Build A22^-1 from pedigree
# ============================================================
cat("\n[3/6] Building A22^-1 from pedigree...\n")
ped <- fread(ped_file, header=FALSE)
setnames(ped, c("id","sire","dam"))
ped[, id   := as.character(id)]
ped[, sire := ifelse(sire=="0", NA, as.character(sire))]
ped[, dam  := ifelse(dam=="0",  NA, as.character(dam))]

ainv_obj <- makeAinv(pedigree = ped)
Ainv_sp  <- ainv_obj$Ainv
Ainv_ids <- rownames(Ainv_sp)

# Match to genotype order
idx <- match(ids, Ainv_ids)
if (any(is.na(idx))) {
  stop("ERROR: Some Ginv IDs not in pedigree A^-1!")
}

A22inv_sp <- Ainv_sp[idx, idx, drop=FALSE]

diagA22inv <- diag(A22inv_sp)
cat("   - mean diag(A22^-1):", mean(diagA22inv), "\n")

# ============================================================
# 4) Scale Ginv to match A22^-1
# ============================================================
cat("\n[4/6] Scaling Ginv to match A22^-1 diagonal mean...\n")
scale_factor <- mean(diagA22inv) / mean(diagG)
cat("   - scale factor =", scale_factor, "\n")

Ginv <- Ginv * scale_factor

diagG2 <- diag(Ginv)
cat("   - mean diag(Ginv) AFTER scaling:", mean(diagG2), "\n")

# ============================================================
# 5) Ensure positive definiteness (fix tiny negatives)
# ============================================================
cat("\n[5/6] Checking for tiny negative diagonals...\n")
neg_idx <- which(diag(Ginv) < 0)
if (length(neg_idx) > 0) {
  cat("   - Found", length(neg_idx), "negative diag entries; setting to small positive value.\n")
  diag(Ginv)[neg_idx] <- 1e-6
} else {
  cat("   - No negative diagonals detected.\n")
}

# Optional: symmetrize (just to be safe)
cat("   - Symmetrizing Ginv...\n")
Ginv <- 0.5 * (Ginv + t(Ginv))
A22inv <- as.matrix(A22inv_sp)   # convert sparse A22^-1 to dense

Ginv_final = 0.95 * Ginv  + 0.05 * A22inv

# ============================================================
# 6) Convert BACK to sparse triplet and export
# ============================================================
cat("\n[6/6] Converting back to sparse and writing output files...\n")
Ginv_sp_new <- Matrix(Ginv_final, sparse = TRUE)
rm(Ginv_final); gc()

Gtrip_new <- summary(Ginv_sp_new)
cat("   - New non-zero entries:", nrow(Gtrip_new), "\n")

# Write .giv
fwrite(Gtrip_new, out_giv, col.names = FALSE, sep=" ")

# Write .order (same IDs)
writeLines(ids, out_order)

cat("\n=====================================================\n")
cat("✅ DONE: Scaled G^-1 saved as", out_giv, "\n")
cat("✅ Order file saved as", out_order, "\n")
cat("=====================================================\n")
cat("Use this in ASReml:\n")
cat("   !WORKSPACE 250000\n")
cat("   !DENSEG\n")
cat("   Ginv_sheep50k_FIXED.giv  !HINV Ginv_sheep50k_FIXED.order\n")
cat("=====================================================\n")
