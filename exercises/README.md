# Exercise/Self-reflecting questions to be completed during the EU-LI-PHE Summer School workshop.

## Part I

### Exercise 1

1. Using the function `get_sensitivity_specificity`, calculate the performance metrics of the IGRA test at the following delta cut off thresholds: -500, -80, -40, 0, 40, 80, 120, 500.

2. Imagine that you are involved in determining and selecting the optimum threshold for the IGRA test within a country. Assuming a disease prevalence of 6%, the potential risks associated with bovine tuberculosis (bTB) disease to humans (e.g., zoonotic), the cost to the national exchequer associated with compensation to farmers who have animals slaughtered that test positive, based on the evidence provided, what would be an optimum threshold to pick for the IGRA test?

3. Assume the prevalence of the disease decreased to 3%, would it be beneficial to increase the stringency of the test or relax it?

### Exercise 2

1. Evaluate the significance of the correlation between _Genotype Principal Component 1_ and _Holstein pedigree percentage_
    - First determine if it should be Pearson or Spearman correlation.
    - perform the calculation with `cor.test`.
    - Plot the values as a scatter plot.

2. How might researchers with only transcriptomics data infer population structure within their data?

### Exercise 3

1. How many genes are significantly differentially expressed at a false discovery rate (FDR) of 0.0001?
    - How many of these genes are exhibiting increased or decreased expression?

2.  The default null hypothesis (Ho) tested in a differential expression analysis is that the logarithmic fold change (log2FC) of a gene between bTB+ and bTB- cattle is exactly 0. Taking guidance from the DESEQ2 manual (https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) Evaluate the hypothesis that the log2FC is > 0.5. Hint: This applied when executing the `results` function of the `dds` object.

3. Many studies conduct a standard DEA as we have done and impose a post log2FC cut-off (e.g., 1, 2 etc.). Is this statistically appropriate to do?

### Exercise 4

1. Perform the `gProfiler` analysis setting the `ordered_query` variable to `F` and determine if we identify any new pathways are identified by conducting an _overrepresentation analysis (ORA)_.

2. Perform the `gProfiler` analysis using only genes exhibiting increased expression and separately decreased expression - you can set the `ordered_query` variable to whatever value you like.
    - Why might it be useful to look at these genes separately?

3. Change the multiple-testing correction method to determine the impact of this method on the number of significant pathways identified.

4.  Perform the `gProfiler` analysis with an input set of more highly significant or less significant DEGs (e.g., padj < 0.00001 | padj < 0.05) and determine if we identify any new pathways are identified

## Part II

### Exercise 1

1. How many _cis_-eQTLs/_cis_-eGenes are identified when correcting for multiple genes tested chromosome wide using;
    - The Bonferroni method
    - The Benjamini-Hochberg method


### Exercise 2

1. What is the distribution of _cis_-eQTLs around the transcriptional start site (TSS) (normal, skewed left, skewed right, bimodal)

2. What is the distribution of the effect size of _cis_-eQTLs (normal, skewed left, skewed right, bimodal)

3. Explain why there are no _cis_-eQTLs were identified with an effect size close to 0

4. What is the median effect size estimate of identified _cis_-eQTLs

5. What is the correlation between effect size and;
    - Distance to TSS
    - Minor allele frequency (MAF)


### Exercise 3

1. Explain how eigen-MT corrects 

2. Why do we need to apply another multiple testing correction across tested genes?


### Exercise 4

1. What is the relationship between the SNP and the expression of _

2. How can this be interpreted in the context of bTB disease

3. Would this SNP be a suitable candidate to breed for/against bTB? Should we assess the functional impact of this using _in vitro_ genome editing?
