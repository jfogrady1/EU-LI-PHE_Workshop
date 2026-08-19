There are two folders in this directory that provide the data for the two practical components

## In the folder `Part_I`, we have the following files

1. `GSE255724_count_matrix_clean.txt` - This matrix contains the expression level of each gene for all _n_ = 123 samples in the dataset. The matrix is organized into an _n x m_ format with _n_ genes and _m_ samples. Each value is an integer, which represents the raw expression of that gene in the sample.

2. `metadata.txt` - This .txt file contains metadata information for all 123 samples including 
- sample ID, 
- Condition [Control (bTB-) | Infected (bTB+)], 
- RNA sequencing batch, 
- Admixture Component 1 & Admixture Component 2 from the `ADMIXTURE` program (https://dalexander.github.io/admixture/) - Note: the sum of both values is equal to 1.
- Eigenvector entries for genotype principal components (PC) 1 and 2 from the `PLINK` program (https://www.cog-genomics.org/plink/)
- Pedigree based Holstein% (i.e., based on familial/historical mating patterns, what % of this animal's genome is purebred holstein. These estimates are not as good as genomic estimates owing to recombination and independent assortment)

3. `IGRA_measurements.txt` - These values represent the change (Δ) in PPD response, calculated by subtracting the optical density (OD) of interferon-γ released by T cells incubated with purified protein derivative avium (PPDa; a cocktail of _Mycobacterium avium_ antigens), measured using an enzyme-linked immunosorbent assay (ELISA), from the OD of interferon-γ released by T cells incubated with purified protein derivative bovine (PPDb; a cocktail of _M. bovis_ antigens).

## In the folder `Part_II`, we have the following files.

1. `ALL.covariates.txt`. These are the covariates used int he mapping of molecular QTLs.
- Transcriptomic principal component 1 (PC1) to 14 - from the `PCA4QTL` R package (https://github.com/heatherjzhou/PCAForQTL). These represent predominantly unknown sources of variation in the data that we need to account for.
- Genotype PC1 to Genotype PC5 from `PLINK`. These represent variance attributed to population/genetic structure that we need to account for.
- Condition
- Age

2. `ALL.expr_tmm_inv.Chr5.bed.gz` - This will be used as input for the molecular QTL mapping with tensorQTL
 - Chromosome
 - start of TSS
 - End of TSS (+ 1bp)
 - gene/phenotype_id (Ensembl)
 - Normalized and standardized expression values for all _n_ = 123 samples
 - The `.gz` extension means that it is compressed. You can uncompress it and looka at the first 10 lines with the following: `gunzip -c ALL.expr_tmm_inv.Chr5.bed.gz | head`
 - We also have the index of this file `ALL.expr_tmm_inv.Chr5.bed.gz.tbi`, which allows programs to find coordinates within the file very fast - these files typically add `.tbi` to the end of the file name.


3. Three `ALL_genotypes_Chr5` files with the extensions `.bed`, `.bim`, and `.fam`, respectively. These are standard binary (`.bed`) representations and ancillary information (`.bim`, and `.fam`) of `PLINK` files. And will be used as input for tensorQTL
 - In many programs, we specify the prefix (i.e., `ALL_genotypes_Chr5`) and the relevant files are collected by whatever program you are using

4. `ALL_genotypes_Chr5.vcf.gz` This file represents a _variant call format *(VCF)*_ file and contains information on all common (MAF > 0.05) variants on BTA 5 for all _n_ = 123 samples
- We also have the index of this file `ALL_genotypes_Chr5.vcf.gz.tbi`, which allows programs to find coordinates within the file very fast - these files typically add `.tbi` to the end of the file name.

5. `Residualised_expression.Chr5.bed.gz`: This is in nenarly the same format as `ALL.expr_tmm_inv.Chr5.bed.gz` but expression values for each gene represent _residuals_ after all the variance from the covariates described in `ALL.covariates.txt` have been removed.
 - We will use these values to plot our _cis_-eQTLs as these values were those that were tested int he association test with tensorqtl.
 - If you are interested, this file was created using the `Qtltools correct` command (https://qtltools.github.io/qtltools/pages/QTLtools-correct.1.html).
 - Again, we also have the index of the file: `Residualised_expression.Chr5.bed.gz.tbi`.