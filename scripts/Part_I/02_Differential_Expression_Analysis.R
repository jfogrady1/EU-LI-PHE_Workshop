##################################################################################
# EU-LI-PHE COST ACTION Summer School Workshop 2026
# Epidemiology, Genetics and Modelling for Infectious Disease Control
# Author: John F. O'Grady
##################################################################################

## Summary - This script will help us identify differentially expressed genes (DEGs) between M. bovis-infected (n = 60) and control non-infected cattle (n = 63)
# It will also help us identify what pathways the significant DEGs are enriched in, which will help us understand the peripheral host response to M. bovis infection in cattle.


####################################################################
# Learning outcomes
####################################################################

# Learn how to conduct a differential expression analysis using DESeq2
# Understand the importance of correcting for population structure in differential expresssion analysis
# Become familiar with multiple hypothesis testing and how to implement appropriate statistical tests
# Learn how to perform gene set enrichment analysis using g:Profiler with appropriate background gene sets
# Become familiar with visualising differential expression results using volcano plots and enrichment results using jitter plots


# Load in the libraries
library(tidyverse)
library(dplyr)
library(DESeq2)
library(ggplot2)
library(apeglm)
library(data.table)
library(gprofiler2)
library(ggrepel)
library(viridis)
library(RColorBrewer)
library(gprofiler2)

# DESeq2 requires 
  # 1. Count matrix (genes are rows, samples are columns)
  # 2. Metadata file  for input into the DESeq2 object - we will use the same file in the population structure analysi

# Read in the data and convert counts to matrix format and metadata to data frame format (required by DESeq2)
counts <- as.matrix(read.table('data/Part_I/GSE255724_count_matrix_clean.txt', skip = 0, header = T, sep="\t",row.names=1)) # matrix file
coldata <- as.data.frame(read.table('data/Part_I/metadata.txt', sep = '\t', skip = 0, header = T, row.names=1)) # metadata file - note th
  
# All rownames of the metadata must match colnames of the count matrix
# This is required by DeSeq2 to ensure that the correct metadata is associated with the correct sample
all(rownames(coldata) == colnames(counts))

# Convert variables to factors and scale continuous variables that will be used in the DESeq2 model.
coldata$Batch <- factor(coldata$Batch, levels = c("1", "2"))
coldata$Condition <- factor(coldata$Condition, levels = c("Control", "Infected"), labels = c("Control", "Infected")) # DESeq2 doesn't like the '+' or '-' sign
coldata$Age <- scale(coldata$Age, center = TRUE)


# Create the DESeq2 object and specify the model
# Note: Variable of interest comes last, in this case Condition
ddsMat <- DESeqDataSetFromMatrix(countData = counts,
                                    colData = coldata,
                                    design = ~ Batch + Geno_PC1 + Geno_PC2 + Age + Condition)

# Filter for low gene counts
# Genes with low counts have high variance and can skew the results (leading to false positives)
# Here we will filter for genes that have at least 6 counts in 20% of the samples
keep <- rowSums(counts(ddsMat) >= 6) >= (ncol(counts) * 0.2) # remove low count genes
ddsMat <- ddsMat[keep,]


###############################################################################
###############################################################################
# Variance stabilising transformation by DESeq2
###############################################################################
###############################################################################


# RNA-seq data is heteroskedastic, meaning that the variance is not constant across the range of mean values. As
# the mean increases so does the variance, which is problematic for downstream analyses that assume homoskedasticity (constant variance). The variance stabilising transformation (VST) 
# is a method that transforms the data to make the variance more constant across the range of mean values. 
# This is important for downstream analyses such as PCA and clustering that are *extremely* sensitive to heteroskedasticity. 

# The mathematics of the model is beyond the scope of this workshop but you can read more about it in the DESeq2 paper and vignette

ddsvsd <- vst(ddsMat)


##################################
##################################
# PCA of the transcriptomics data
##################################
##################################


# If we conduct a PCA on all genes, many will be correlated,
# We only want to visualise the top n genes, this is arbitrary but 1,500 are often used in the literature. 
#This is a good compromise between including enough genes to capture the biological signal and not including too many genes that may introduce noise.

# perform PCA with the plotPCA function from DESeq2
pcaData <- plotPCA(ddsvsd, intgroup=c("Condition"), ntop = 1500, returnData=TRUE) # top 1500 variable genes

# Store the proportion of variation explained by each PC in a varible
percentVar <- round(100 * attr(pcaData, "percentVar"))

# We will show the impact of population structure by layering over other metadata available to us
# This information comes from our population structure analysis
pcaData$Holstein_admixture <- coldata$Admixture_Component_1 



# Here add that bins Admixture proportions into suitable bins for visualisation. This is arbitrary but we will use 0-1% (none), 1-20%, 20-40%, 40-60%, 60-80% and 80-100% bins.
pcaData <- pcaData %>% mutate(ADMIX_bin = cut((Holstein_admixture * 100), right = F, breaks=c(0.000,1,20, 40, 60, 80, 100), 
    labels = c("0", "1-20", "20-40", "40-60", "60-80", "80-100")))
colnames(pcaData)[7] <- "%Holstein"


# We will plot eigenvector 1 (PC1) and eigenvector 2 (PC2) and colour by condition and shape by admixture proportion
# This will be saved as variable 'pca'
pca <- ggplot(pcaData, aes(PC1, PC2, color=Condition, shape = `%Holstein`)) +
geom_point(size=4, alpha = 0.8) +
scale_color_manual(values = c("#2166ac", "#b2182b")) +
scale_shape_manual(values = c(15,16,17,18,11,13)) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
labs(shape = 'Admixture\ncomponent 1') + theme_classic(base_size = 14) +
  # Make it look better with the (theme() function)
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 14),
    legend.text = element_text(size = 14),
    legend.position = "right"
  )

pca

###############################################################################
###############################################################################
# Differential Expression Anlaysis
###############################################################################
###############################################################################


# Perform the differential expression
dds <- DESeq(ddsMat)


# Extract the results from the DESeq2 object and perform log fold change shrinkage using the apeglm method. 
# This helps to reduce the noise in the log fold change estimates, especially for genes with low counts or high variance. 
res <- results(dds, contrast=c("Condition","Control","Infected"))
res <- lfcShrink(dds, coef="Condition_Infected_vs_Control", type="apeglm")

# Get a quick overview of the results
summary(res, alpha = 0.05)
str(res)

# Save as a data frame
res_df <- as.data.frame(res)

res_df$padj <- as.numeric(res_df$padj)    

min(res_df$log2FoldChange, na.rm = TRUE)
# Get the gene symbols for the Ensembl IDs
# Here we are using the gconvert function from gprofiler
all_symbols <- gconvert(query = rownames(res_df), organism = "btaurus", 
        target="ENSG", mthreshold = Inf, filter_na = FALSE)


# If the symbol does not exist, replace with the default ENSGID
# This can happen if ensembl has updated ENSGID names
res_df$Symbol <- all_symbols$name 
res_df$Symbol[is.na(res_df$Symbol)] <- rownames(res_df)[is.na(res_df$Symbol)] # replace NA with the Ensembl ID


# Will create a variable to colour our points
res_df$diffexpressed <- "Not DE"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
res_df$diffexpressed[res_df$log2FoldChange > 0 & res_df$padj < 0.05] <- "DE Up"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
res_df$diffexpressed[res_df$log2FoldChange < 0 & res_df$padj < 0.05] <- "DE Down"

# Make a variable to plot some key genes
res_df <- res_df %>% mutate(PLOT_Symbol =
                            case_when(
                                log2FoldChange < -.1 & padj < 0.00005 ~ Symbol,
                                log2FoldChange > 1 & padj < 0.00005 ~ Symbol,
                                log2FoldChange > 0 & padj < 0.0000005 ~ Symbol,
                                log2FoldChange < -0.5 & padj < 0.01 ~ Symbol,
                                log2FoldChange > 1.5 & padj < 0.001 ~ Symbol,
                                padj < 0.00005 ~ Symbol,
                                FALSE ~ ""))

# table counts the number of entries for each category in a vector
# We will save the results of DE up and DE down for annotation on the volcano plot
tab <- table(res_df$diffexpressed)
DE_down <- tab[1]
DE_up <- tab[2]

# If the gene symbol starts with ENSBTAG, it is an uncharacterised gene 
# and we will not plot it on the volcano plot.
res_df <- res_df %>%
  mutate(PLOT_Symbol = if_else(
    grepl("^ENSBTAG", PLOT_Symbol),
    NA_character_,
    PLOT_Symbol
  ))


# Here we will create a volcano plot to visualise the results
Volcano <- ggplot(data=res_df, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed, label=PLOT_Symbol)) +
geom_point(size = 1.5, alpha = 0.5) + # increased point size and added alpha for transparency
scale_color_manual("Comparison", values=c("#2166ac", "#b2182b", "grey")) +
labs(x=expression(log[2]("fold change")),
    y=expression(-log[10](italic(P)[adj])))+
scale_x_continuous(limits = c(-4,4), breaks = c(-4,-3,-2,-1,0,1,2,3,4)) +
scale_y_continuous(limits = c(0,8.2), breaks = c(0,1,2,3,4,5,6,7,8)) +
geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed") +
geom_text_repel(colour = "black", fontface = 4, max.overlaps = 40, size = 3.5) +
scale_y_continuous(expand = c(0, 0), limits = c(0, 8.5)) +
#scale_x_continuous(expand = c(0, 0), limits = c(-1, 1)) +
theme_classic(base_size = 14) +
  # Make it look better with the (theme() function)
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 14),
    legend.text = element_text(size = 14),
    legend.position = "right"
  ) +
annotate("text", x=-3, y=-log10(1e-8), size = 5, label=paste(DE_down), col="#2166ac", fontface = "bold", alpha = 0.8) +
annotate("text", x=-2.66, y=-log10(1e-8), label= "↓", col="#2166ac", size = 10, fontface = "bold", alpha = 0.8) +
annotate("text", x=3, y=-log10(1e-8), size = 5, label=paste(DE_up), col="#b2182b", fontface = "bold", alpha = 0.8) +
annotate("text", x=2.66, y=-log10(1e-8), label="\u2191", col="#b2182b", size = 10, fontface = "bold", alpha = 0.8) +
guides(color = guide_legend(override.aes = list(size = 4)))

Volcano


res_final <- res_df %>%
filter(padj <= 0.05)
res_final <- res_final %>% select(-c("diffexpressed", "PLOT_Symbol"))
res_df <- res_df %>% filter(!is.na(padj)) # Only include genes tested
dim(res_df)
res_final <- res_final[order(res_final$padj, decreasing = F),]
dim(res_final)

head(res_final)
head(res_df)




#############################################
# Exercises
#1. How many genes are significantly differentially expressed at a false discovery rate (FDR) of 0.0001?
    #- How many of these genes are exhibiting increased or decreased expression?
    # Hint: go back to the res_final data frame and filter for padj < 0.0001 and then count the number of genes with log2FoldChange > 0 and log2FoldChange < 0.

#2. The default null hypothesis (Ho) tested in a differential expression analysis is that the logarithmic fold change (log2FC) of a gene between bTB+ and bTB- cattle is exactly 0. Taking guidance from the DESEQ2 manual (https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) Evaluate the hypothesis that the log2FC is > 0.5. Hint: This applied when executing the `results` function of the `dds` object.

#3. Many studies conduct a standard DEA as we have done and impose a post log2FC cut-off (e.g., 1, 2 etc.). Is this statisticaly appropriate to do?

#############################################
















#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################

###############################################################################
###############################################################################
# Functional Enrichment Analysis of DEGs using g:Profiler
###############################################################################
###############################################################################


# we will reduce the number of genes to highly significant DEGs (padj < 0.01) for functional enrichment analysis
# This will make the computational time slightly faster and 'should' improve the quality of the results.
res_filtered <- res_final %>% filter(padj < 0.01)
genes <- rownames(res_filtered)

# convert to a list
# Gprofiler accepts this format as input
genes <- list(genes)

# gost is a function from the gprofiler2 package that performs functional enrichment analysis using the g:Profiler web service.
# Results will change from time to time as the g:Profiler database is updated regularly.
# We are performing an enrichment based on P-value
results <- gost(query = genes,organism = "btaurus", correction_method = "fdr", ordered_query = T, domain_scope = "known", custom_bg = rownames(res_df), user_threshold = 0.05, sources = c("GO:BP", "GO:CC", "KEGG", "REAC"), evcodes = T)

# Here is a vector of some key terms we might be interested in highlighting on the plot
# This is arbitrary and can be changed to whatever terms you wish to highlight
terms <- c("MDA-5 signaling pathway",
        "defense response to virus",
        "Interferon Signaling",
        "cytokine production",
        "type I interferon-mediated signaling pathway",
        "Influenza A",
        "RAF3-dependent IRF activation pathway",
        "interleukin-27-mediated signaling pathway",
        "interferon-mediated signaling pathway", "innate immune response",
        "MyD88-independent toll-like receptor signaling pathway",
        "regulation of cytokine-mediated signaling pathway",
        "type II interferon-mediated signaling pathway",
        "programmed necrotic cell death",
        "NOD-like receptor signaling pathway",
        "RIG-I-like receptor signaling pathway",
        "ISG15 antiviral mechanism",
        "Yersinia infection",
        "regulation of transformation of host cell by virus",
        "Antiviral mechanism by IFN-stimulated genes",
        "ISG15 antiviral mechanism ",
        "complex of collagen trimers",
        "NF-kB activation through FADD/RIP-1 pathway mediated by caspase-8 and -10",
        "Host-pathogen interaction of human coronaviruses - interferon induction",
        "Toll Like Receptor 4 (TLR4) Cascade",
        "Coronavirus disease - COVID-19",
        "IRF3-mediated induction of type I IFN",
        "Integrin cell surface interactions"
        )

# Term names are in the column result under the results gprofiler object
# Here we are creating a column to label if we want to label it.
results$result <- results$result %>%
mutate(label = ifelse(term_name %in% terms, term_name, NA))


#####################################################
#####################################################
# Jitter plot of functional enrichment results
#####################################################
#####################################################

# Make a new dataframe with the results of interest
result_df <- results$result

# Create an alpha value to highlight our terms of interest
result_df <- result_df %>% mutate(alpha_value = ifelse(term_name %in% terms, 1, 0.25))

# Colour palette from the viridis R package
# There are lots of colour palettes out there
# We will use the package RcolourBrewer to select a colour palette.
# Usually, I like to use the following website to get an idea of suitable colours: https://colorbrewer2.org/#type=sequential&scheme=BuGn&n=3

# Keep the position of the points consistent for reproducibility
pos <- position_jitter(width = 0.3, seed = 3)

# Now make the plot
my_enrichment_plot <- ggplot(result_df, aes(x = source, y = -log10(p_value))) +
geom_jitter(aes(color = source), alpha = result_df$alpha_value, position = pos, size = 3) +  # Use shape instead of color
scale_color_brewer(palette = "Dark2" , name = "Source") +  # Apply the color palette
scale_y_continuous(limits = c(0, 10), breaks = seq(0:10)) +
labs(y = expression(-log[10](italic(P)[adj])), x = "") +
theme_classic(base_size = 14) +
  # Make it look better with the (theme() function)
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 14),
    legend.text = element_text(size = 14),
    legend.position = "right"
  ) +
geom_text_repel(aes(x=source, y=-log10(p_value)), label = results$result$label, max.overlaps = 20,
                size=4.0, color='black', fontface = "bold",
                fill='#FFFFFF33',
                position = pos,
                force = 4
) +
geom_hline(yintercept=-log10(0.05), col="black", linetype = "dashed") +
guides(color = guide_legend(override.aes = list(size = 5)))


library(patchwork)
Volcano + my_enrichment_plot + plot_annotation(tag_levels = 'A')


##########################################################################################
##########################################################################################
# 1. Perform the GProfiler analysis setting the `ordered_query` variable to `F` and determine if we identify any new pathways are identified by conducting an _overrepresentation analysis (ORA)_.

# 2. Perform the Gprofiler analysis with an input set of more highly significant or less significant DEGs (e.g., padj < 0.00001 OR padj < 0.05) and determine if we identify any new pathways are identified

# 3. Perform the Gprofiler analysis using only genes exhibiting increased expression and seperately decreased expression - you can set the `ordered_query` variable to whatever value you like.
   # - Why might it be useful to look at these genes seperately?

# 4. Change the multiple-testing correction method to determine the impact of this method on the number of significant pathways identified.
   # - You can look at the options in the gprofiler2 R package documentation (https://cran.r-project.org/web/packages/gprofiler2/gprofiler2.pdf) for more information on the different methods available.
##########################################################################################
##########################################################################################

