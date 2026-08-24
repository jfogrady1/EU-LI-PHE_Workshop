## Load in libraries
library(data.table)
library(tidyverse)
library(ggplot2)


# Read in the data
# fread comes from the data.table package
IFNG_measurements <- fread("data/Part_I/IGRA_measurements.txt", header = TRUE, sep = "\t")



# Plot the data in ggplot as a boxplot and separate by the condition (Control vs Infected) and add a horizontal line at 80 (the threshold for a positive result)
# Note the value of 80 (purple line) corresponds to the cut-off used for this assay by the Department of Agriculture Food and the Marine.
ggplot(aes(x = Condition, y = delta_PPD), data = IFNG_measurements) +
  geom_boxplot(aes(fill = Condition), alpha = 0.5, outlier.colour = NA) +
  labs(x = "Condition", y = "Delta PPD") +
  # Here we are overlaying the points as a jitter (not just a straight line)
  # Setting the position of the points with seed so that they don't change and move around
   geom_jitter(aes(color = Condition), position = position_jitter(seed = 12, width = 0.2), size = 3) +
  theme_bw(base_size = 14, base_family = "Helvetica") +
  scale_fill_manual(values = c("Control" = "#1f78b4", "Infected" = "#b2182b")) +
  scale_color_manual(values = c("Control" = "#1f78b4", "Infected" = "#b2182b")) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "right"
  ) + geom_hline(yintercept = 80, linetype = "dashed", color = "purple", linewidth = 1.5)



# Look at the distribution of the data
# Can visually assess if it is normally distributed
# This will help determine what statistical test we will use to compare the two groups
# Pretty obvious though that the distribution of the two groups is significantly different, with the control group being skewed to the left and the infected group being skewed to the right.
ggplot(data = IFNG_measurements, aes(x = delta_PPD, fill = Condition)) + facet_wrap(~Condition) + geom_histogram(alpha = 0.5) +
  theme_classic(base_size = 14, base_family = "Helvetica") +
  scale_fill_manual(values = c("Control" = "#1f78b4", "Infected" = "#b2182b")) +
   theme(
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text = element_text(colour = "black", size = 12),
    axis.text.x = element_text(colour = "black", size = 12, angle = 45, hjust = 1, vjust = 1),
    axis.title = element_text(colour = "black", size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "top"
  )

# Here we will determine if the data is normally distributed using the Shapiro-Wilk test (shapiro.test() function)
# Note: the '%>%' operator lets us 'pipe' the output from one function into another
IFNG_measurements %>% group_by(Condition) %>% 
summarise(shapiro_p = shapiro.test(delta_PPD)$p.value, 
            mean = mean(delta_PPD), 
            sd = sd(delta_PPD), 
            median = median(delta_PPD), 
            IQR = IQR(delta_PPD))



# Can see that they are not normally distributed (particularly in  control group), so we will use a non-parametric test (Wilcoxon rank sum test) to compare the two groups
wilcox.test(delta_PPD ~ Condition, data = IFNG_measurements)

# Plot the value again
ggplot(aes(x = Condition, y = delta_PPD), data = IFNG_measurements) +
  geom_boxplot(aes(fill = Condition), alpha = 0.5, outlier.colour = NA) +
  labs(x = "Condition", y = "Delta PPD") +
  geom_jitter(aes(color = Condition), position = position_jitter(seed = 12, width = 0.2), size = 3) +
  theme_bw(base_size = 14, base_family = "Helvetica") +
  scale_fill_manual(values = c("Control" = "#1f78b4", "Infected" = "#b2182b")) +
  scale_color_manual(values = c("Control" = "#1f78b4", "Infected" = "#b2182b")) +
  theme(
    axis.line = element_line(colour = "black", linewidth = 0.6),
    axis.ticks = element_line(colour = "black", linewidth = 0.6),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(colour = "black", size = 14),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    legend.position = "right"
  ) + geom_hline(yintercept = 80, linetype = "solid", color = "purple", linewidth = 1.5, linewidth = 0.6)



# Assess the sensitivity and specificity of the assay using the get_sensitivity_specificity() function
# This is the structure for a function in R
get_sensitivity_specificity <- function(values,
                                        labels,
                                        threshold = 80,
                                        positive_label = Infected,
                                        na.rm = FALSE) {

  positive_label <- 'Infected'
  negative_label <- 'Control'

  pred_positive <- values >= threshold

  actual_positive <- labels == positive_label
  actual_negative <- labels == negative_label

  tp <- sum(pred_positive & actual_positive)
  fn <- sum(!pred_positive & actual_positive)
  tn <- sum(!pred_positive & actual_negative)
  fp <- sum(pred_positive & actual_negative)

  sensitivity <- tp / (tp + fn)
  specificity <- tn / (tn + fp)

  # Produce the outputs as a list.
  list(
    sensitivity = sensitivity,
    specificity = specificity,
    confusion = c(TP = tp, FN = fn, TN = tn, FP = fp),
    threshold = threshold,
    positive_label = positive_label,
    n = length(values)
  )
}

# Apply the function to our data.
get_sensitivity_specificity(values = IFNG_measurements$delta_PPD, 
labels = IFNG_measurements$Condition, 
threshold = 80, 
positive_label = "Infected", 
na.rm = TRUE)  
