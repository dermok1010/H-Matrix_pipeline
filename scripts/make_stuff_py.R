#!/usr/bin/env Rscript
# ============================================================
# export_for_python.R  (clean, lossless, no dense CSVs)
#
# Prepares inputs for the Python H^-1 builder:
#   - Builds A^-1 (full pedigree) with nadiv (dgCMatrix)
#   - Extracts A22^-1 for genotyped animals (same sparse type)
#   - Exports both in Matrix Market (.mtx) to preserve structure
#   - Writes order files for A and G blocks
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(nadiv)
})

# -------------------------------
# 1) Arguments
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript export_for_python.R <pedigree.csv> <geno_id_file> <output_prefix>")
}

ped_file      <- args[1]
geno_id_file  <- args[2]
output_prefix <- args[3]

cat("\n=================================================\n")
cat("  export_for_python.R\n")
cat("  Pedigree:    ", ped_file, "\n")
cat("  Geno IDs:    ", geno_id_file, "\n")
cat("  Output pref: ", output_prefix, "\n")
cat("=================================================\n")

# -------------------------------
# 2) Read pedigree
# -------------------------------
cat("\n[1/4] Reading pedigree...\n")
ped <- fread(ped_file, header = FALSE)
if (ncol(ped) < 3) stop("Pedigree must have 3 columns: id sire dam")
setnames(ped, c("id","sire","dam"))
ped[, id   := as.character(id)]
ped[, sire := fifelse(sire=="0", NA_character_, as.character(sire))]
ped[, dam  := fifelse(dam=="0",  NA_character_, as.character(dam))]
cat("   - Pedigree rows:", nrow(ped), "\n")

# -------------------------------
# 3) Read genotyped IDs
# -------------------------------
cat("\n[2/4] Reading genotyped IDs...\n")
geno_ids <- fread(geno_id_file, header = FALSE)[[1]]
geno_ids <- trimws(as.character(geno_ids))
cat("   - Genotyped animals:", length(geno_ids), "\n")

# -------------------------------
# 4) Build A^-1 (full pedigree)
# -------------------------------
cat("\n[3/4] Building A^-1 (sparse dgCMatrix) via nadiv...\n")
# [3/4] Build A^-1 (sparse) and force order to match the input pedigree
cat("\n[3/4] Building A^-1 (sparse dgCMatrix) via nadiv...\n")
ainv_obj <- makeAinv(pedigree = ped)
Ainv_sp  <- ainv_obj$Ainv                       # dgCMatrix
Ainv_ids <- rownames(Ainv_sp)

# --- force Ainv to EXACTLY match input pedigree order (ped$id) ---
# ped$id is in .SRT order if you passed the SRT file
perm <- match(ped$id, Ainv_ids)
if (any(is.na(perm))) {
  miss <- Ainv_ids[is.na(match(Ainv_ids, ped$id))]
  stop(paste("Internal error: Ainv had IDs not found in input pedigree order. First few:",
             paste(head(miss), collapse=", ")))
}
Ainv_sp <- Ainv_sp[perm, perm, drop = FALSE]
Ainv_ids <- ped$id  # now Ainv ids == input order
dimnames(Ainv_sp) <- list(Ainv_ids, Ainv_ids)

nA <- length(Ainv_ids)
cat("   - A^-1 dimension:", nA, "x", nA, " | nnz:", length(Ainv_sp@x), "\n")
cat("   - Order enforced to match input pedigree (first 5): ",
    paste(head(Ainv_ids, 5), collapse=", "), "\n")

# Extract A22^-1 using G IDs IN THIS ORDER
idxG <- match(geno_ids, Ainv_ids)
if (any(is.na(idxG))) {
  missing <- geno_ids[is.na(idxG)]
  stop(paste("ERROR: Some genotyped IDs not in pedigree (first few):",
             paste(head(missing), collapse=", ")))
}
A22inv_sp <- Ainv_sp[idxG, idxG, drop = FALSE]
cat("   - A22^-1 dimension:", paste(dim(A22inv_sp), collapse=" x "),
    " | nnz:", length(A22inv_sp@x), "\n")

if (any(is.na(idxG))) {
  missing <- geno_ids[is.na(idxG)]
  stop(paste("ERROR: Some genotyped IDs not in pedigree (first few):",
             paste(head(missing), collapse=", ")))
}
A22inv_sp <- Ainv_sp[idxG, idxG, drop = FALSE]
cat("   - A22^-1 dimension:", paste(dim(A22inv_sp), collapse=" x "),
    " | nnz:", length(A22inv_sp@x), "\n")

# -------------------------------
# 5) Export files for Python
# -------------------------------
cat("\n[4/4] Exporting for Python (Matrix Market)...\n")

# (a) Full A^-1: Matrix Market
ainv_mtx <- paste0(output_prefix, "_Ainv.mtx")
writeMM(Ainv_sp, ainv_mtx)

# (b) A22^-1: Matrix Market
a22_mtx <- paste0(output_prefix, "_A22inv.mtx")
writeMM(A22inv_sp, a22_mtx)

# (c) Order files
writeLines(Ainv_ids, paste0(output_prefix, "_Ainv.order"))
writeLines(geno_ids, paste0(output_prefix, "_G_ids.txt"))

cat("\n✅ Export complete\n")
cat("   - ", ainv_mtx, "\n", sep = "")
cat("   - ", a22_mtx, "\n", sep = "")
cat("   - ", output_prefix, "_Ainv.order\n", sep = "")
cat("   - ", output_prefix, "_G_ids.txt\n", sep = "")
cat("=================================================\n")
