
# Load the libraries
library(data.table)
library(tidyverse)
library(ggplot2)


# Read in the eQTL data
eQTL_data <- fread("results/Part_II/ALL.cis_qtl_fdr0.05.txt")

# Look at the structure
str(eQTL_data)
head(eQTL_data, 10)

# How many eGenes
table(eQTL_data$is_eGene)

# Filter for significant eQTLs
eQTL_data_sig <- eQTL_data %>% filter(is_eGene == TRUE)

######################################
# Properties of cis-eQTLs.
######################################

# Here we will plot the distribution of the distance to the TSS for the top cis-eQTLs for each eGene.
ggplot(eQTL_data_sig, aes(x = start_distance)) +
  geom_histogram(alpha = 0.5, bins = 30, col = "#322525", fill = 'steelblue') +
  labs(x = "Distance to Transcriptional Start Site (TSS)", y = "Count") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 14)
  )

# Questions to consider:
# Where do cis-eQTLs tend to be located relative to the TSS of their associated gene? 
# Are they enriched in promoters, enhancers, or other regulatory regions?


# Here we will plot the distribution of the effect sizes (slopes) for the top cis-eQTLs for each eGene.
ggplot(eQTL_data_sig, aes(x = slope)) +
  geom_histogram(alpha = 0.5, bins = 30, col = "#322525", fill = 'steelblue') +
  labs(x = "Effect Size (Slope)", y = "Count") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 14)
  )

# Questions to consider:
# What is the distribution of effect sizes for the top cis-eQTLs?
# Why are there no eQTLs with effect sizes close to 0?




# Here we will plot the distribution of the absolute effect sizes (slopes) for the top cis-eQTLs for each eGene.
ggplot(eQTL_data_sig, aes(y = abs(slope), x = 'Top cis-eQTLs')) +
  geom_boxplot(alpha = 0.5, col = "#322525", fill = 'steelblue', outlier.colour = NA) +
  theme_classic(base_size = 14) +
  labs (x = "Top cis-eQTLs", y = "Absolute effect size (slope)") +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 14)
  ) + geom_jitter(position = position_jitter(seed = 12, width = 0.2), size = 3)


summary(abs(eQTL_data_sig$slope))

# Questions to consider:
# What is the distribution of absolute effect sizes for the top cis-eQTLs?
# What is a limitation of our data?



# Here we will plot the relationship between the absolute effect sizes (slopes) for the top cis-eQTLs for each eGene
# and the distance to the TSS
ggplot(eQTL_data_sig, aes(x = abs(start_distance), y = abs(slope))) +
  geom_point(alpha = 0.5, col = "#322525", fill = 'steelblue') +
  labs(x = "Distance to Transcriptional Start Site (TSS)", y = "Absolute effect size (slope)") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 14)
  ) + geom_smooth(method = "lm", col = "red", se = FALSE)

# Questions to consider:
# What is the relationship between effect size and distance to TSS? Are larger effect sizes associated with variants closer to the TSS?



# Here we want to see the relationship between absolute effect size and the minor allele frequency (MAF) of the top cis-eQTLs for each eGene.
# Here we will calculate the MAF.
# Why are we using 246?
eQTL_data_sig$MAF <- eQTL_data_sig$ma_count / 246


## Do the plotting.
ggplot(eQTL_data_sig, aes(y = MAF, x = abs(slope))) +
  geom_point(alpha = 0.5, col = "#322525", fill = 'steelblue') +
  labs(x = "Absolute effect size (slope)", y = "Minor Allele Frequency (MAF)") +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 14)
  ) + geom_smooth(method = "lm", col = "red", se = FALSE) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black", linewidth = 0.6)




# Let's get a sense of the number of variants tested in each cis window.
ggplot(eQTL_data_sig, aes(y = num_var, x = 'Group')) +
  geom_boxplot(alpha = 0.5, col = "#322525", fill = 'steelblue', outlier.colour = NA) +
  labs(y = "Number of variants in cis-window", x = "Count") +
  theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 14)
  ) + geom_jitter(position = position_jitter(seed = 12, width = 0.2), size = 3)

summary(eQTL_data_sig$num_var)



# Determine if there are any pleiotropic eQTLs.

# Filter for variants found in multiple genes
pleio <- eQTL_data_sig %>%
  group_by(variant_id) %>%
  summarise(
    n_genes = n_distinct(phenotype_id),
    genes   = paste(phenotype_id, collapse = ","),
    .groups = "drop"
  ) %>%  arrange(desc(n_genes))


# Plot the results as a barplot for each gene.
ggplot(data.frame(pleio), aes(x = n_genes)) +
  geom_bar(fill = "steelblue") +
  scale_x_continuous(breaks = 1:max(pleio$n_genes)) +
  labs(
    x = "Number of eGenes sharing the same lead variant",
    y = "Number of variants",
    title = "Distribution of cis-eQTL pleiotropy"
  ) + theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 14)
  )


######################################
# Visualization of cis-eQTLs.
######################################

## Plot some eQTLs of interest

# To plot eQTL data we need genotype information (from vcf)
# We also need residualized expression information (.txt file)


# vcf data
library(vcfR)
vcf_data <- read.vcfR("data/Part_II/ALL_genotypes_Chr5.vcf.gz")


vcf_geno <- extract.gt(vcf_data, element = "GT", as.numeric = FALSE)
colnames(vcf_geno) <- gsub(".*_","", colnames(vcf_geno))

vcf_geno <- gsub("0/0", "0", vcf_geno)
vcf_geno <- gsub("0/1", "1", vcf_geno)
vcf_geno <- gsub("1/0", "1", vcf_geno)
vcf_geno <- gsub("1/1", "2", vcf_geno)

# Expression information 
residualized_expr <- fread("data/Part_II/Residualised_expression.Chr5.bed.gz")
residualized_expr <- residualized_expr %>% select(-c("#chr", "start", "end", "id", "strd")) %>% column_to_rownames(var = "gid")

residualized_expr <- residualized_expr %>% filter(rownames(residualized_expr) %in% eQTL_data_sig$phenotype_id)

all_symbols <- gconvert(query = rownames(residualized_expr), organism = "btaurus", 
        target="ENSG", mthreshold = Inf, filter_na = FALSE)
head(all_symbols)
all(rownames(residualized_expr) == all_symbols$input)

all_symbols <- all_symbols %>% mutate(name = case_when(is.na(name) ~ input, TRUE ~ name)) %>% select(input, name)

rownames(residualized_expr) <- all_symbols$name
eQTL_data_sig <- eQTL_data_sig %>% mutate(phenotype_id = case_when(phenotype_id %in% all_symbols$input ~ all_symbols$name, TRUE ~ phenotype_id))
# Function
plot_genotype_boxplot <- function(expr_mat, geno_mat, gene, variant,
                                   geno_labels = c("0", "1", "2"),
                                   log_expr = FALSE,
                                   jitter = TRUE,
                                   point_alpha = 0.6,
                                   colors = c("#4E79A7", "#F28E2B", "#E15759")) {
 
  if (!gene %in% rownames(expr_mat))
    stop(sprintf("Gene '%s' not found in rownames(expr_mat)", gene))
  if (!variant %in% rownames(geno_mat))
    stop(sprintf("Variant '%s' not found in rownames(geno_mat)", variant))
 
  # match samples between the two matrices
  common_samples <- intersect(colnames(expr_mat), colnames(geno_mat))
 
  expr_vals <- as.numeric(expr_mat[gene, common_samples])
  geno_vals <- as.numeric(geno_mat[variant, common_samples])
 
  df <- data.frame(
    sample     = common_samples,
    expression = expr_vals,
    genotype   = factor(geno_vals, levels = c(0, 1, 2), labels = geno_labels)
  )
  df <- df[!is.na(df$genotype) & !is.na(df$expression), ]
 
  n_labels <- table(df$genotype)
  x_labels <- paste0(levels(df$genotype), "\n(n=", n_labels, ")")
 
  p <- ggplot(df, aes(x = genotype, y = expression, fill = genotype)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6, width = 0.6) +
    labs(
      x = paste0("Genotype of variant ", variant),
      y = paste0("Residualized Expression of gene ", gene)
    ) +
    scale_x_discrete(labels = x_labels) +
    scale_fill_manual(values = colors) +
    theme_classic(base_size = 14) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 14)
  )
 
  if (jitter) {
    p <- p + geom_jitter(width = 0.15, size = 1.5, alpha = point_alpha)
  }
 
  p
}

head(eQTL_data_sig)
head(residualized_expr)

# Filter for gene of interest
eQTL_data_sig[eQTL_data_sig$phenotype_id == 'LRRK2', ]



# Recent evidence suggests that LRRK2 is a negative regulator of M. tuberculosis phagosome maturation in macrophages.
# We have identified a cis-eQTL associated with the expression of this gene, and inhibiting LRRK2 activity in mice increased the ability
# of macrophages to control M. tuberculosis growth.
# Let's see what we can visualize with this.
#https://pmc.ncbi.nlm.nih.gov/articles/PMC6003659/

plot_genotype_boxplot(expr_mat = residualized_expr, 
                      geno_mat = vcf_geno, 
                      gene = "LRRK2", 
                      variant = "5:40592048:A:G")


head(eQTL_data_sig[eQTL_data_sig$phenotype_id == 'CLEC4E', ])

# CLEC4E is a C-type lectin receptor that recognizes mycobacterial components and is involved in the innate immune response to M. tuberculosis infection.
# It works in combination with TLR4 and is required for the induction of autophagy to restrict and destroy M. tuberculosis in macrophages.
#(https://www.tandfonline.com/doi/full/10.1080/15548627.2019.1658436#abstract)

# We have identified a cis-eQTL for this gene.
plot_genotype_boxplot(expr_mat = residualized_expr, 
                      geno_mat = vcf_geno, 
                      gene = "CLEC4E", 
                      variant = "5:101690540:T:C")


eQTL_data_sig

# Let's look at a negative association
plot_genotype_boxplot(expr_mat = residualized_expr, 
                      geno_mat = vcf_geno, 
                      gene = "TBC1D15", 
                      variant = "5:1444559:G:T")
