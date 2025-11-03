#!/usr/bin/env Rscript
# ============================================================
# make_Ginv_ASReml.R (FINAL, CLEAN)
# Build G^-1 (.giv + .order) for ASReml using:
#  - pedigree (A22)
#  - genomic relationship (G)
# Proper pipeline:
#   1) A^-1 from pedigree (nadiv)
#   2) Extract A22^-1, invert to A22
#   3) Build raw G from PLINK
#   4) match.G2A to align A22 & G
#   5) G.tuneup align + blend(0.95/0.05) + bend
#   6) G.inverse(sparse)
#   7) Export .giv + .order
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(nadiv)
  library(ASRgenomics)
})

# -------------------------------
# 1) Input arguments
# -------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript make_Ginv_ASReml.R <pedigree_csv> <GRM_prefix> <output_prefix>")
}

ped_file      <- args[1]
grm_prefix    <- args[2]
output_prefix <- args[3]

blend_weight <- 0.95  # 95% G, 5% A22

cat("\n=================================================\n")
cat("  make_Ginv_ASReml.R (FINAL)\n")
cat("  Pedigree:    ", ped_file, "\n")
cat("  GRM prefix:  ", grm_prefix, "\n")
cat("  Output pref: ", output_prefix, "\n")
cat("  Blend:       ", blend_weight, "G +", 1-blend_weight, "A\n")
cat("=================================================\n")

# -------------------------------
# 2) Read pedigree (id, sire, dam), unknown=0
# -------------------------------
cat("\n[1/7] Reading pedigree...\n")
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
cat("\n[2/7] Reading genotyped IDs...\n")
geno_ids <- fread(paste0(grm_prefix, ".grm.id"), header = FALSE)[[1]]
geno_ids <- as.character(geno_ids)
geno_ids <- trimws(geno_ids)
nG <- length(geno_ids)
cat("   - Genotyped animals:", nG, "\n")

# -------------------------------
# 4) Build A^-1 sparse via nadiv, extract A22^-1
# -------------------------------
cat("\n[3/7] Building full A^-1 (sparse) with makeAinv...\n")
ainv_obj <- makeAinv(pedigree = ped)
Ainv_sp  <- ainv_obj$Ainv  # sparse inverse matrix
Ainv_ids <- rownames(Ainv_sp)

# Match genotyped IDs to A^-1 order
idxG <- match(geno_ids, Ainv_ids)
if (any(is.na(idxG))) {
  stop("ERROR: Some genotyped IDs not in pedigree A^-1. Fix pedigree or IDs.")
}
cat("   - A^-1 dimension:", dim(Ainv_sp), "\n")
cat("   - Genotyped IDs found in A^-1:", sum(!is.na(idxG)), "\n")

# Extract A22^-1
A22inv_sp <- Ainv_sp[idxG, idxG, drop=FALSE]
cat("   - A22^-1 dimension:", dim(A22inv_sp), "\n")

# -------------------------------
# 5) Invert A22^-1 -> full A22
# -------------------------------
cat("\n[4/7] Inverting A22^-1 to full A22 (dense)...\n")
A22_full <- solve(A22inv_sp)  # ~13GB matrix
rm(A22inv_sp); gc()

# Assign names
rownames(A22_full) <- geno_ids
colnames(A22_full) <- geno_ids

# Ensure base 'matrix' class
A22_full <- as.matrix(A22_full)
cat("   - A22 dimension:", dim(A22_full), "\n")

# -------------------------------
# 6) Build raw G (full) from PLINK
# -------------------------------
cat("\n[5/7] Building raw G from PLINK .grm...\n")
Gtab <- fread(paste0(grm_prefix, ".grm"), header = FALSE)
setnames(Gtab, c("i","j","rel","Nsnps"))

Gmat <- matrix(0, nG, nG)
for (r in seq_len(nrow(Gtab))) {
  i <- Gtab$i[r]; j <- Gtab$j[r]
  Gmat[i,j] <- Gtab$rel[r]
  Gmat[j,i] <- Gtab$rel[r]
}

rownames(Gmat) <- geno_ids
colnames(Gmat) <- geno_ids
Gmat <- as.matrix(Gmat)
cat("   - Gmat dimension:", dim(Gmat), "\n")

# -------------------------------
# 7) match.G2A to align G & A22
# -------------------------------
cat("\n[6/7] Matching G with A22 via match.G2A()...\n")
matched <- match.G2A(
  A = A22_full,
  G = Gmat,
  clean = TRUE,
  ord = TRUE,
  mism = TRUE
)
A22clean <- matched$Aclean
Gclean   <- matched$Gclean

stopifnot(all(rownames(A22clean) == rownames(Gclean)))
cat("   - Post-match dimension:", dim(Gclean), "\n")

# -------------------------------
# 8) Tune G: align -> blend -> bend
# -------------------------------
cat("\n[7/7] Tuning G (align -> blend -> bend)...\n")

# 1) ALIGN
cat("   - Align G to A22...\n")
G1 <- G.tuneup(
  G = Gclean,
  A = A22clean,
  align = TRUE,
  blend = FALSE,
  bend = FALSE
)$Gb

# 2) BLEND (0.95 G + 0.05 A22)
cat("   - Blend 95% G + 5% A22...\n")
G2 <- G.tuneup(
  G = G1,
  A = A22clean,
  align = FALSE,
  blend = TRUE,
  pblend = 1 - blend_weight, # 0.05
  bend = FALSE
)$Gb

# 3) BEND (ensure positive-definite)
cat("   - Bend G to ensure PD...\n")
Gfinal <- G.tuneup(
  G = G2,
  bend = TRUE,
  eig.tol = 1e-6,
  align = FALSE,
  blend = FALSE
)$Gb

cat("   - Gfinal diag: min =", min(diag(Gfinal)),
    " mean =", mean(diag(Gfinal)),
    " max =", max(diag(Gfinal)), "\n")

# -------------------------------
# 9) Invert tuned G -> sparse G^-1
# -------------------------------
cat("\n[8/8] Inverting tuned G to sparse G^-1...\n")
Ginv_out <- G.inverse(Gfinal, sparseform = TRUE)
Ginv_trip <- Ginv_out$Ginv
cat("   - Ginv non-zeros:", nrow(Ginv_trip), "\n")

# -------------------------------
# 10) Export .giv + .order
# -------------------------------
giv_file   <- paste0(output_prefix, ".giv")
order_file <- paste0(output_prefix, ".order")

fwrite(Ginv_trip, giv_file, col.names = FALSE, sep = " ")
writeLines(rownames(Gfinal), order_file)

cat("\n=================================================\n")
cat("✅ DONE: G^-1 written to ", giv_file, "\n", sep="")
cat("✅ Order file written to ", order_file, "\n", sep="")
cat("Use in ASReml:\n")
cat("  !GIV ", basename(giv_file), " ", basename(order_file), "\n", sep="")
cat("=================================================\n")
