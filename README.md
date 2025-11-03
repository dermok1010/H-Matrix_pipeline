# H-Matrix_pipeline
This pipeline is for constructing, and validating G and H inverse matrices for mixed-model genetic evaluation in ASReml

## Scripts

build_ginv.R
Constructs a genomic inverse relationship matrix that is aligned, blended, and positive definite.
The script reads a pedigree and PLINK-formatted relationship matrix, and outputs an ASReml-ready .giv and .order files

fix_Ginv.R
This script ensures the diagnonal mean and scale of the G inverse matrix matches that of the A22 inverse (the pedigree-based relationship for genotyped animals).
The script reads the .giv and .order files that were created in build_ginv.R and also the pedigree file

lets_make_Hinv.py.py
Constructs a combined pedigree and genomic inverse matrix (H inverse matrix)
Requires:
Sparse full pedigree A inverse matrix
Sparse A22 inverse
.giv file
.order file

Outputs:
....._Hinv.giv Unbent H inverse
....._Hinv.order full pedigree ID used in H inverse
....._Hinv_PD.giv Positive definite (ridge-bent) version of H inverse
....._Hinv_PD.order matching order for PD matrix

## Slurm Scripts
Example slurm scripts are shown in slurm_scripts folder, these are guidelines and should be modified to suit your own system/paths
