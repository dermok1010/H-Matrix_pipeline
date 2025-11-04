# H-Matrix_pipeline
This pipeline is for constructing, and validating G and H inverse matrices for mixed-model genetic evaluation in ASReml

## Scripts

### build_ginv.R
Constructs the genomic inverse relationship matrix (**G⁻¹**) that is aligned with the pedigree, blended for scale stability, and positive definite.

It reads a **sorted three-column pedigree file** (`id`, `sire`, `dam`) and a **PLINK-formatted genomic relationship matrix** (`.grm` and `.grm.id`).

The pedigree **must be sorted**, meaning parents must appear before offspring (as in ASReml `!TRIM !SORT` format).  
This ensures correct construction of A⁻¹ and A₂₂⁻¹ and alignment with the genomic data.

The script aligns G with A₂₂, applies a blend (typically 95% G and 5% A₂₂), bends the result to ensure positive-definiteness, and outputs two ASReml-ready files:  
a sparse genomic inverse (`.giv`) and a matching ID order file (`.order`).

---

### fix_Ginv.R
Ensures the diagonal scale and mean of the genomic inverse match those of the pedigree-based A₂₂ inverse.

It reads the `.giv` and `.order` files produced by `build_ginv.R`, along with the same pedigree file.  
The script rebuilds A₂₂⁻¹ from the pedigree, rescales and symmetrises the G⁻¹ matrix as needed, and writes a corrected version of the genomic inverse with the same sparse `.giv` / `.order` format.

---

### make_stuff_py.R
Prepares pedigree information for the Python stage of the pipeline.

It reads the **sorted pedigree file** and the `.grm.id` list of genotyped animals, then builds sparse full-pedigree (**A⁻¹**) and genotyped-subset (**A₂₂⁻¹**) inverses using `nadiv`.  
Both are exported in **Matrix Market (.mtx)** format, along with two text files: one listing the full pedigree order and one listing the genotyped IDs.  
These provide all pedigree components required for building H⁻¹.

---

### lets_make_Hinv_py.py
Constructs the combined pedigree and genomic inverse (**H⁻¹**) used in single-step genomic evaluation.

It requires five inputs: the full A⁻¹ `.mtx`, the A₂₂⁻¹ `.mtx`, the genomic inverse `.giv`, its `.order` file, and an output prefix.  
The script embeds the genomic correction term (G⁻¹ − A₂₂⁻¹) into A⁻¹ to form H⁻¹, checks symmetry and scale, and if necessary adds a small ridge to guarantee positive-definiteness.  
It outputs both an unbent and a bent (PD) version of the H inverse, each with corresponding `.giv` and `.order` files ready for ASReml.

---

### diagnostic_hinv.py
Performs structural diagnostics on any `.giv` file—typically the H⁻¹ produced in the previous step.

It checks matrix symmetry, diagonal range, sparsity, and positive-definiteness, reporting a concise summary to the terminal.  
This script is mainly for verification, ensuring that the matrix is well-behaved and suitable for use in ASReml or downstream mixed-model analyses.
## Slurm Scripts
Example slurm scripts are shown in slurm_scripts folder, these are guidelines and should be modified to suit your own system/paths
