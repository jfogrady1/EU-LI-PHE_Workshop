There are two folders in this directory that provide the data for the two practical components

In the folder `Part_I`, we have the following files

1. `GSE255724_count_matrix_clean.txt` - This matrix contains the expression level of each gene for all _n_ = 123 samples in the dataset. The matrix is organized into an _n x m_ format with _n_ genes and _m_ samples. Each value is an integer, which represents the raw expression of that gene in the sample.

2. `metadata.txt` - This .txt file contains metadata information for all 123 samples including 
- sample ID, 
- Condition [Control (bTB-) | Infected (bTB+)], 
- RNA sequencing batch, 
- Admixture Component 1 & Admixture Component 2 from the `ADMIXTURE` program (https://dalexander.github.io/admixture/) - Note: the sum of both values is equal to 1.
- Eigenvector entries for genotype principal components (PC) 1 and 2 from the `PLINK` program (https://www.cog-genomics.org/plink/)
- Pedigree based Holstein% (i.e., based on familial/historical mating patterns, what % of this animal's genome is purebred holstein. These estimates are not as good as genomic estimates owing to recombination and independent assortment)

3. `IGRA_measurements.txt` - These values represent the change (Δ) in PPD response, calculated by subtracting the optical density (OD) of interferon-γ released by T cells incubated with purified protein derivative avium (PPDa; a cocktail of _Mycobacterium avium_ antigens), measured using an enzyme-linked immunosorbent assay (ELISA), from the OD of interferon-γ released by T cells incubated with purified protein derivative bovine (PPDb; a cocktail of _M. bovis_ antigens).
