# Load in the libraries
library(data.table)
library(tidyverse)
library(ggplot2)
library(gprofiler2)
library(ggpubr)
library(ggdist)
library(gghalves)

# Filter for significant ieQTLs
ieQTL_results <- fread("results/Part_II/ALL_TB_interaction.cis_qtl_top_assoc.txt.gz")


# Inspect the data
str(ieQTL_results)
head(ieQTL_results, 10)

# Count the number of significant ieQTLs
ieQTL_sig <- ieQTL_results %>% filter(pval_adj_bh < 0.05) %>% nrow()

# Save to a new data frame
ieQTL_sig <- ieQTL_results %>% filter(pval_adj_bh < 0.05)
ieQTL_sig

# Get the gene name
all_symbols <- gconvert(query = ieQTL_sig$phenotype_id, organism = "btaurus", 
        target="ENSG", mthreshold = Inf, filter_na = FALSE)

ieQTL_sig$name <- all_symbols$name


# Load the VCF data
library(vcfR)
vcf_data <- read.vcfR("data/Part_II/ALL_genotypes_Chr5.vcf.gz")

vcf_geno <- extract.gt(vcf_data, element = "GT", as.numeric = FALSE)
colnames(vcf_geno) <- gsub(".*_","", colnames(vcf_geno))
vcf_geno <- gsub("0/0", "0", vcf_geno)
vcf_geno <- gsub("0/1", "1", vcf_geno)
vcf_geno <- gsub("1/0", "1", vcf_geno)
vcf_geno <- gsub("1/1", "2", vcf_geno)
# Load the expression information again
residualized_expr <- fread("data/Part_II/Residualised_expression.Chr5.bed.gz")
residualized_expr <- residualized_expr %>% select(-c("#chr", "start", "end", "id", "strd")) %>% column_to_rownames(var = "gid")
residualized_expr <- residualized_expr %>% filter(rownames(residualized_expr) %in% ieQTL_sig$phenotype_id)
all(rownames(residualized_expr) == all_symbols$input)
rownames(residualized_expr) <- all_symbols$name
ieQTL_sig <- ieQTL_sig %>% mutate(phenotype_id = case_when(phenotype_id %in% all_symbols$input ~ all_symbols$name, TRUE ~ phenotype_id))
head(ieQTL_sig, 10)
common_samples <- intersect(colnames(residualized_expr), colnames(vcf_geno))

# Create the data frame
df <- data.frame(
    sample = common_samples,
    expression = as.numeric(residualized_expr["PCBP2", common_samples]),
    genotype = factor(vcf_geno["5:25861761:G:A", common_samples],
                      levels = c(0,1,2)),
    condition = c(rep('Control', 63), rep('Infected', 60))
)


# Plot the data
ggplot(df, aes(y = expression, x = genotype, shape = condition, colour = condition)) + 
    ggdist::stat_halfeye(data = subset(df, genotype == '0' & condition =="Control"), aes(x = genotype, y = expression), slab_colour =  "#2166ac", slab_fill = "#2166ac", adjust = .33, width = .33, alpha=0.5, position = position_nudge(x=0), side = "left") +
    ggdist::stat_halfeye(data = subset(df, genotype == '0' & condition =="Infected"), aes(x = genotype, y = expression), slab_colour =  "#b2182b", slab_fill = "#b2182b", adjust = .33, width = .33, alpha=0.5, position = position_nudge(x=0), side = "right") +
    ggdist::stat_halfeye(data = subset(df, genotype == '1' & condition =="Control"), aes(x = genotype, y = expression), slab_colour =  "#2166ac", slab_fill = "#2166ac", adjust = .33, width = .33, alpha=0.5, position = position_nudge(x=0), side = "left") +
    ggdist::stat_halfeye(data = subset(df, genotype == '1' & condition =="Infected"), aes(x = genotype, y = expression), slab_colour = "#b2182b", slab_fill = "#b2182b", adjust = .33, width = .33, alpha=0.5, position = position_nudge(x=0), side = "right") +
    ggdist::stat_halfeye(data = subset(df, genotype == '2' & condition =="Control"), aes(x = genotype, y = expression), slab_colour =  "#2166ac", slab_fill = "#2166ac", adjust = .33, width = .33, alpha=0.5, position = position_nudge(x=0), side = "left") +
    ggdist::stat_halfeye(data = subset(df, genotype == '2' & condition =="Infected"), aes(x = genotype, y = expression), slab_colour = "#b2182b", slab_fill = "#b2182b", adjust = .33, width = .33, alpha=0.5, position = position_nudge(x=0), side = "right") +
    xlab('5:25861761:G:A') + ylab(paste0("Residualised expression of ", "PCBP2")) + labs(fill = 'Genotype') +
    scale_x_discrete(labels =  c('0', '1', '2')) +
    geom_point(size = 4) +
    scale_colour_manual(values = c("#2166ac","#b2182b")) +
    geom_smooth(method = "lm", aes(group = condition, colour = condition), se = FALSE) +
    theme(axis.text.x = element_text(angle = 0, size = 15, colour = "black"),
    axis.text.y = element_text(angle = 0, size = 15, colour = "black"),
    axis.title.y = element_text(size = 21, color = "black"),
    axis.title.x = element_text(size = 21, color = "black"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(size = 15, color = "black"),
    legend.text = element_text(size = 15)) +  theme_classic(base_size = 14, base_family = "Helvetica") +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 14),
    legend.text = element_blank(),
    legend.title = element_blank(),
    legend.position = "none"
  )

# PCBP2 binds to poly(C) stretches of DNA and RNA and negatively regulates cyclic GMP-AMP synthase (cGAS) activity.
# Cyclic-GAS is a cytosolic DNA sensor that, upon binding to DNA, converts ATP and GTP to cyclic GMP-AMP (cGAMP), which subsequently binds to and activates STING, thereby triggering type I interferon production.
# For tuberculous mycobacteria, cGAS and subsequent intracellular innate immune response mechanisms are activated by the mycobacterial ESX-1 secretion system.
# Overexpression of PCBP2 has been shown to reduce cGAS-STING signalling, whereas the opposite pattern has been observed following PCBP2 loss.


# What does this mean for bTB?