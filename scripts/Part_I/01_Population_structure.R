## Summary - This script will help us get an understanding of the genetic structure present in
# our bovine tuberculosis (bTB) disease dataset of M. bovis-infected (n = 60) and control non-infected cattle (n = 63) 
# that we will be investigating as part of this workshop. This is important as inter- and 
# intra-breed genomic variation exists, which can impact the host response to infectious disease challenge,
# can confound genome-wide association studies and as such, needs to be controlled for and appreciated.

# Load in libraries
library(tidyverse)
library(data.table)
library(ggplot2)

# Read in the data
data <- read.table('data/Part_I/metadata.txt', header = TRUE, sep = '\t')


# Inspect the data
str(data)

# Look at first 10-rows
head(data, 10)

# Make sure all samples were loaded in
dim(data)


# Some data manipulation
data$Condition <- factor(data$Condition, levels = c("Control", "Infected"), labels = c('bTB-', 'bTB+'))



#########
# Investigating population structure
#########


# PCA plot of genotype data
# PCA generated from ~23,000 SNPs using PLINK v1.9 (https://www.cog-genomics.org/plink/1.9/) 
# and the --pca function, which uses the EIGENSTRAT algorithm to perform PCA on genotype data.
p1 <- ggplot(data, aes(x = Geno_PC1, y = Geno_PC2, color = Condition)) +
  geom_point(size = 3) +
  labs(x = "Principal component 1 (PC1)",
       y = "Principal component 2 (PC2)") +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_classic(base_size = 14) +
  scale_color_manual(values = c("bTB-" = "#1f78b4", "bTB+" = "#b2182b")) +
  # Make it look better with the (theme() function)
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "right"
  )

p1

# Structure plot of genotype data
# Admixture values generated using ADMIXTURE v1.3.0 (https://dalexander.github.io/admixture/) and the --cv flag to determine the optimal K-value (K = 2).

# Need to format the data into long format (currently in wide format)
# So we can plot the two components as a stacked bar plot

data_admixture <- data %>%
  pivot_longer(
    cols = starts_with("Admixture_Component"),
    names_to = "K",
    values_to = "Proportion"
  )


sample_order <- data_admixture %>%
  filter(K == "Admixture_Component_1") %>%
  arrange(Condition, desc(Proportion)) %>% 
  pull(Sample)

data_admixture <- data_admixture %>%
  mutate(Sample = factor(Sample, levels = sample_order))

  
p2 <- ggplot(data_admixture,
       aes(x = Sample,
           y = Proportion,
           fill = K)) +
  geom_col(width = 1) +
  labs(x = "Sample ID",
       y = "Admixture proportion (K = 2)") +
  scale_fill_manual(values = c(
    "Admixture_Component_1" = "#364C54",
    "Admixture_Component_2" = "#94475E"
  )) + theme_classic(base_size = 14) +
  scale_y_continuous(expand = c(0, 0)) +
  # Make it look better with the (theme() function)
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text = element_text(colour = "black", size = 12),
    axis.text.x = element_text(colour = "black", size = 12, angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_text(colour = "black", size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "top"
  ) + geom_vline(xintercept = 63.5, linetype = "dashed", color = "black", size = 0.6, linewidth = 2)

p2

# Do ADMIXTURE and PCA display report compatible information

# We will investigate the relationship between the admixture component 1 and PC1
# The easiest way to assess the relationship between variables x and y is via correlation.


# lets look at the distribution of the two variables
# Assessing the distribution of the data is the first task in any analysis.
# If the data is normally distributed (bell shaped) we can apply parametric tests (e.g., Pearson correlation), if not, we can apply non-parametric tests (e.g., Spearman correlation) or alternatively
# transform the data to a normal distribution e.g., via log transformation

# We will do this in base R as these plots don't need to be fancy
hist(data$Admixture_Component_1)
hist(data$Geno_PC1)

# Dont appear to be normally distributed but can foramally test with shapiro.wilk test
shapiro.test(data$Admixture_Component_1)
shapiro.test(data$Geno_PC1)

# A significant p-value (p < 0.05) indicates that the data is not normally distributed, which is the case for both variables here. 
# As such, we will use a non-parametric test (Spearman) to assess the correlation between the two variables.


admix_pc_cor <- cor.test(data$Admixture_Component_1, data$Geno_PC1, method = "spearman", exact = FALSE)
admix_pc_cor

# report exact p-value
admix_pc_cor$p.value

# This is a highly significant correlation, lets see how it looks visually
p3 <- ggplot(data, aes(x = Admixture_Component_1, y = Geno_PC1)) +
  geom_point(col = 'darkgrey', size = 3) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.6) +
  labs(x = "Admixture component 1",
       y = "Principal component 1 (PC1)") +
  theme_classic(base_size = 14) +
  # Make it look better with the (theme() function)
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "right"
  ) + annotate('text', x = 0.15, y = 0.1, label = paste0("Spearman correlation = ", round(admix_pc_cor$estimate, 3)), size = 3)
p3

# PCA and ADMIXTURE information are (relatively) mutually compatible, at least for this ancestry component (K = 2). This is important as it means that the two methods are capturing the same underlying genetic structure in our dataset, which is reassuring.
# However, what is this ancestry component? is it a breed, geographic location, familial, sire-based?
# We hypothesise that principal component of genomic variation (i.e., PC1) can be explained by Holstein breed percentage

hist(data$Holstein_percentage)
shapiro.test(data$Holstein_percentage)


admix_holstein_cor <- cor.test(data$Admixture_Component_1, data$Holstein_percentage, method = "spearman", exact = FALSE)
admix_holstein_cor
admix_holstein_cor$p.value

# lets plot this
p4 <- ggplot(data, aes(y = Holstein_percentage, x = Admixture_Component_1)) +
  geom_point(col = 'darkgrey', size = 3) +
  labs(y = "Holstein breed percentage",
       x = "Admixture component 1")  +
  theme_classic(base_size = 14) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  # Make it look better with the (theme() function)
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "right"
  ) + geom_smooth(method = "lm", se = TRUE, color = "steelblue", linewidth = 1) +
  geom_abline(slope = 100, intercept = 0, linetype = "dashed", color = "darkred", size = 0.6, linewidth = 2) +
  annotate('text', x = 0.15, y = 85, label = paste0("Spearman correlation = ", round(admix_holstein_cor$estimate, 3)), size = 3)
p4



# Lets make a publication ready figure of the population structure analyses
library(patchwork)

# Combine the plots into a single figure
combined_plot <- (p1 + p3 + p4) / p2 + plot_annotation(tag_levels = 'A')
combined_plot

# Need to change around the tags
# Because we want admixture plot to be panel B
p1 <- p1 + labs(tag = "A")
p2 <- p2 + labs(tag = "B")
p3 <- p3 + labs(tag = "C")
p4 <- p4 + labs(tag = "D")

combined_plot <- (p1 + p3 + p4) / p2
combined_plot


# save the plot
# Will make it extra wide for the admixture plot
# If you want to save this, uncomment the line below but create the directory results/figures/ first
#ggsave("results/figures/population_structure.pdf", combined_plot, width = 20, height = 10, units = "in", dpi = 600)

## evaluate the correlation between genoPC1 and Holstein%.

