install.packages("depmixS4")

library(depmixS4)
library(tidyverse)

# Calculate average expression for all genes
average_expression <- AverageExpression(rgc, assay = 'RNA', group.by = 't_factor')

# The result is a list, with one component for each assay. Extract the RNA assay component.
avg_expr_matrix <- average_expression$RNA

timepoint_split <- SplitObject(rgc, split.by = 't_factor')

# Initialize a list to store standard deviations
std_expr_list <- list()

# Loop through each timepoint and calculate standard deviation for each gene
for (timepoint in names(timepoint_split)) {
  # Extract the RNA data for the current timepoint
  expr_data <- GetAssayData(timepoint_split[[timepoint]], assay = "RNA", slot = "data")
  
  # Calculate the standard deviation for each gene (row) across cells (columns)
  std_expr_list[[timepoint]] <- apply(expr_data, 1, sd)
}

# Combine the results into a matrix
std_expr_matrix <- do.call(cbind, std_expr_list)
colnames(std_expr_matrix) <- names(std_expr_list)

# Filter out genes where all expression values are zero across all timepoints
non_zero_genes <- rowSums(avg_expr_matrix != 0) > 0

# Subset both the average expression matrix and the standard deviation matrix
avg_expr_matrix_filtered <- avg_expr_matrix[non_zero_genes, ]
std_expr_matrix_filtered <- std_expr_matrix[non_zero_genes, ]


# Find the timepoint with the highest average expression for each gene in the filtered matrix
max_timepoints_filtered <- apply(avg_expr_matrix_filtered, 1, function(x) names(which.max(x)))

# Create a data frame with the results for filtered genes
max_expr_df_filtered <- data.frame(
  Gene = rownames(avg_expr_matrix_filtered),
  MaxTimepoint = max_timepoints_filtered,
  MaxExpression = apply(avg_expr_matrix_filtered, 1, max),
  StdDev = mapply(function(gene, timepoint) std_expr_matrix_filtered[gene, timepoint], 
                  rownames(avg_expr_matrix_filtered), max_timepoints_filtered)
)
gc()

# Perform pairwise t-tests for each gene in the filtered data, comparing the max timepoint with others
# Perform pairwise t-tests for each gene in the filtered data, comparing the max timepoint with others
# Function to determine significance based on the new criteria
max_expr_df_filtered$significant <- mapply(function(gene) {
  
  # Get the expression values and standard deviations for the gene across timepoints
  expr_values <- avg_expr_matrix_filtered[gene, ]
  std_devs <- std_expr_matrix_filtered[gene, ]
  
  # Sort the timepoints by expression value (from highest to lowest)
  sorted_timepoints <- names(sort(expr_values, decreasing = TRUE))
  
  # Get the highest and second-highest expression values
  max_timepoint <- sorted_timepoints[1]
  second_max_timepoint <- sorted_timepoints[2]
  
  max_expr <- expr_values[max_timepoint]
  second_max_expr <- expr_values[second_max_timepoint]
  
  max_std_dev <- std_devs[max_timepoint]
  second_max_std_dev <- std_devs[second_max_timepoint]
  
  # Check if the highest value minus its std dev is greater than the second highest value plus its std dev
  if ((max_expr - max_std_dev) > (second_max_expr + second_max_std_dev)) {
    return(TRUE)  # Significant
  } else {
    return(FALSE)  # Not significant
  }
  
}, rownames(avg_expr_matrix_filtered))

# Check the results
head(max_expr_df_filtered)

# Function to determine significance based on the new criteria without lower or upper bound considerations
max_expr_df_filtered$significant_nolowerorupper <- mapply(function(gene) {
  
  # Get the expression values and standard deviations for the gene across timepoints
  expr_values <- avg_expr_matrix_filtered[gene, ]
  std_devs <- std_expr_matrix_filtered[gene, ]
  
  # Sort the timepoints by expression value (from highest to lowest)
  sorted_timepoints <- names(sort(expr_values, decreasing = TRUE))
  
  # Get the highest and second-highest expression values
  max_timepoint <- sorted_timepoints[1]
  second_max_timepoint <- sorted_timepoints[2]
  
  max_expr <- expr_values[max_timepoint]
  second_max_expr <- expr_values[second_max_timepoint]
  
  max_std_dev <- std_devs[max_timepoint]
  second_max_std_dev <- std_devs[second_max_timepoint]
  
  # Check if either condition is satisfied
  if ((max_expr > second_max_expr + second_max_std_dev) || (max_expr - max_std_dev > second_max_expr)) {
    return(TRUE)  # Significant
  } else {
    return(FALSE)  # Not significant
  }
  
}, rownames(avg_expr_matrix_filtered))

# Check the results
head(max_expr_df_filtered)

# Define the list of genes
genes_to_filter <- c("GDF9", "GDF10", "GDF11", "GDF15", "WNT5A", "WNT5B", "WNT11", "NRG1", "NRG2", "NRG3", "NRG4", 
                     "FGF1", "FGF2", "FGF4", "FGF5", "FGF6", "FGF7", "FGF3", "FGF10", "FGF22", "FGF8", "FGF17", 
                     "FGF18", "FGF9", "FGF20", "FGF16", "FGF19", "FGF21", "FGF23", "VEGFA", "VEGFB", "VEGFC", "VEGFD", 
                     "PGF", "DHH", "IHH", "SHH", "MIF", "PTN", "PRSS3", "CTSG", "F2", "PLG", "GZMA", "PRSS2", "PRSS1", 
                     "ADCYAP1", "ENHO", "LAMA1", "LAMA2", "LAMA3", "LAMA4", "LAMA5", "LAMB1", "LAMB2", "LAMB3", 
                     "LAMC1", "LAMC2", "LAMC3", "TNC", "TNR", "TNN", "TNXB", "ALCAM", "APP", "CADM1", "CADM3", 
                     "NFASC", "CNTN2", "CNTN1", "EFNA1", "EFNA2", "EFNA3", "EFNA4", "EFNA5", "EFNB1", "EFNB2", "EFNB3", 
                     "NEGR1", "LRRC4C", "LRRC4", "LRRC4B", "DLK1", "DLL1", "DLL3", "DLL4", "JAG1", "JAG2", "NRXN1", 
                     "NRXN2", "NRXN3", "OCLN", "SEMA7A", "PTPRS", "PPIA", "PCDHA1", "PCDHA2", "PCDHB1", "PCDHB10", 
                     "PCDHGA1", "PCDHGB1", "FLRT1", "FLRT2", "FLRT3", "LRFN4", "LRFN5", "SLC32A1", "GAD1", "SLC6A1", 
                     "GAD2", "SHMT1", "SLC6A5", "SLC6A9", "SHMT2", "GLS", "SLC17A6", "SLC1A1", "CBLN1", "PGD2", 
                     "PTGDS", "PTGES", "PTGES2", "PGE2", "AKR1C3", "PGF2a", "TXA2", "CBR1", "PRXL2B", "CEL", 
                     "DHCR24", "DHCR7", "DHT", "SRD5A1", "AKR1D1", "HSD17B12", "HSD17B3")

# Filter the dataframe by the list of genes
filtered_df <- max_expr_df_filtered[max_expr_df_filtered$Gene %in% genes_to_filter, ]
ordered_df <- filtered_df[match(genes_to_filter, filtered_df$Gene), ]

# Display the ordered dataframe
print(ordered_df)


table(max_expr_df_filtered$MaxTimepoint)
# Visualization with filtered data
gene_of_interest <- "SULT2B"
plot_data_filtered <- data.frame(
  Timepoint = factor(colnames(avg_expr_matrix_filtered), levels = c('t1]','t2]','t3]','t4]','t5]','t6]','t7')),#c('t1]','t2]','t3]','t4]','t5]','t6]','t7')), correct_order
  AverageExpression = avg_expr_matrix_filtered[gene_of_interest, ],
  StdDev = std_expr_matrix_filtered[gene_of_interest, ]
)

# Plot the average expression with error bars
ggplot(plot_data_filtered, aes(x = Timepoint, y = AverageExpression)) +
  geom_line(group = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = AverageExpression - StdDev, ymax = AverageExpression + StdDev), width = 0.2) +
  labs(title = paste("Expression of", gene_of_interest, "Across Timepoints"),
       x = "Timepoint",
       y = "Average Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


all_genes_freq <- as.data.frame(table(max_expr_df_filtered$MaxTimepoint))
significant_genes_freq <- as.data.frame(table(max_expr_df_filtered$MaxTimepoint[max_expr_df_filtered$significant == TRUE]))
significant_nolowerorupper_genes_freq <- as.data.frame(table(max_expr_df_filtered$MaxTimepoint[max_expr_df_filtered$significant_nolowerorupper == TRUE]))

# Rename columns for clarity
colnames(all_genes_freq) <- c("Timepoint", "Frequency")
colnames(significant_genes_freq) <- c("Timepoint", "Frequency")
colnames(significant_nolowerorupper_genes_freq) <- c("Timepoint", "Frequency")

# Convert Timepoint to a factor with the specified order
all_genes_freq$Timepoint <- factor(all_genes_freq$Timepoint, levels = correct_order) #correct order or t1, t2

# Create the plot with the ordered x-axis
ggplot(all_genes_freq, aes(x = Timepoint, y = Frequency)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  labs(title = "Frequency of Max Expression Timepoints (All Genes)",
       x = "Timepoint", y = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

significant_genes_freq$Timepoint <- factor(significant_genes_freq$Timepoint, levels =  correct_order)

# Plot 2: Frequency of MaxTimepoint where significant == TRUE
ggplot(significant_genes_freq, aes(x = Timepoint, y = Frequency)) +
  geom_bar(stat = "identity", fill = "lightgreen") +
  labs(title = "Frequency of Max Expression Timepoints (Significant Genes)",
       x = "Timepoint", y = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

significant_nolowerorupper_genes_freq$Timepoint <- factor(significant_nolowerorupper_genes_freq$Timepoint, levels = correct_order)

# Plot 3: Frequency of MaxTimepoint where significant_nolowerorupper == TRUE
ggplot(significant_nolowerorupper_genes_freq, aes(x = Timepoint, y = Frequency)) +
  geom_bar(stat = "identity", fill = "lightcoral") +
  labs(title = "Frequency of Max Expression Timepoints (Significant No Lower/Upper Genes)",
       x = "Timepoint", y = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

correct_order <- c("Day 59", "Week 9", "Week 11", "Day 80", "Day 82", "Week 12", 
                   "Week 13", "Week 14", "Week 15","Day 105", "Week 16",  
                   "Week 17","Day 125", "Week 18",  "Week 19", "Week 20", 
                   "Week 22", "Week 24", "Week 27")

# Reorder the columns of the average expression matrix
avg_expr_matrix <- avg_expr_matrix[, correct_order]
dim(avg_expr_matrix)
# This matrix now has genes as rows and timepoints as columns
# Remove genes where expression is zero across all timepoints
avg_expr_matrix_filtered <- avg_expr_matrix[rowSums(avg_expr_matrix != 0) > 0, ]
# Filter out genes with zero variance
avg_expr_matrix_filtered <- avg_expr_matrix_filtered[apply(avg_expr_matrix_filtered, 1, var) > 0, ]
# Filter out genes that have fewer than, say, 3 non-zero values across all timepoints
min_non_zero_values <- 7
avg_expr_matrix_filtered <- avg_expr_matrix_filtered[rowSums(avg_expr_matrix_filtered != 0) >= min_non_zero_values, ]

threshold <- 1e-5
avg_expr_matrix_filtered <- avg_expr_matrix_filtered[apply(avg_expr_matrix_filtered, 1, var) > threshold, ]

# Check the dimensions after filtering
dim(avg_expr_matrix_filtered)

markers %>%
group_by(cluster) %>%
slice_max(n=25, order_by = avg_log2FC)



average_expression <- AverageExpression(rgc, assay = 'RNA', group.by = 't_factor') #switch to timepoint for next iteration
avg_expr_matrix <- average_expression$RNA
non_zero_genes <- rowSums(avg_expr_matrix != 0) > 0

# Subset both the average expression matrix and the standard deviation matrix
avg_expr_matrix_filtered <- avg_expr_matrix[non_zero_genes, ]
avg_expr_matrix_filtered




results <- data.frame(
  gene = rownames(avg_expr_matrix_filtered),
  multiple_r_squared = numeric(nrow(avg_expr_matrix_filtered)),
  adjusted_r_squared = numeric(nrow(avg_expr_matrix_filtered)),
  p_value = numeric(nrow(avg_expr_matrix_filtered))
)

# Loop through each gene and fit the model
for (i in 1:nrow(avg_expr_matrix_filtered)) {
  gene_data <- data.frame(
    expression = avg_expr_matrix_filtered[i, ],
    index = 1:ncol(avg_expr_matrix_filtered)
  )
  
  # Fit the polynomial model
  poly_fit <- lm(expression ~ poly(index, 2), data = gene_data)
  
  # Get the summary of the model
  model_summary <- summary(poly_fit)
  
  # Extract multiple R-squared, adjusted R-squared, and p-value
  multiple_r_squared <- model_summary$r.squared
  adjusted_r_squared <- model_summary$adj.r.squared
  p_value <- pf(model_summary$fstatistic[1], model_summary$fstatistic[2], model_summary$fstatistic[3], lower.tail = FALSE)
  
  # Store the results
  results$multiple_r_squared[i] <- multiple_r_squared
  results$adjusted_r_squared[i] <- adjusted_r_squared
  results$p_value[i] <- p_value
}
results

results$significant <- ifelse(results$p_value < 0.05, "significant", "not significant")

table(results$significant)

results_t_factor <- results

correct_order <- c("Day 59", "Week 9", "Week 11", "Day 80", "Day 82", "Week 12", 
                   "Week 13", "Week 14", "Week 15","Day 105", "Week 16",  
                   "Week 17","Day 125", "Week 18",  "Week 19", "Week 20", 
                   "Week 22", "Week 24", "Week 27")

# Reorder the columns of the average expression matrix to match the correct order of timepoints
avg_expr_matrix_filtered <- avg_expr_matrix_filtered[, correct_order]

# Initialize the results data frame
results <- data.frame(
  gene = rownames(avg_expr_matrix_filtered),
  multiple_r_squared = numeric(nrow(avg_expr_matrix_filtered)),
  adjusted_r_squared = numeric(nrow(avg_expr_matrix_filtered)),
  p_value = numeric(nrow(avg_expr_matrix_filtered)),
  stringsAsFactors = FALSE
)

# Loop through each gene and fit the model
for (i in 1:nrow(avg_expr_matrix_filtered)) {
  gene_data <- data.frame(
    expression = avg_expr_matrix_filtered[i, ],
    index = 1:ncol(avg_expr_matrix_filtered)  # Treat the correct_order as the index
  )
  
  # Fit the polynomial model
  poly_fit <- lm(expression ~ poly(index, 2), data = gene_data)
  
  # Get the summary of the model
  model_summary <- summary(poly_fit)
  
  # Extract multiple R-squared, adjusted R-squared, and p-value
  multiple_r_squared <- model_summary$r.squared
  adjusted_r_squared <- model_summary$adj.r.squared
  p_value <- pf(model_summary$fstatistic[1], model_summary$fstatistic[2], model_summary$fstatistic[3], lower.tail = FALSE)
  
  # Store the results
  results$multiple_r_squared[i] <- multiple_r_squared
  results$adjusted_r_squared[i] <- adjusted_r_squared
  results$p_value[i] <- p_value
}

# Add the 'significant' column based on the p-value
results$significant <- ifelse(results$p_value < 0.05, "significant", "not significant")

# View the results
results

table(results$significant)

# Load ggplot2 for plotting
library(ggplot2)

# Filter results for significant genes
significant_results <- results_t_factor[results_t_factor$significant == "significant", ] #or regular results for timepoint

# Extract the adjusted R-squared values
adjusted_r_squared_values <- results$adjusted_r_squared #or significant_results

# Create the plot
ggplot(data.frame(adjusted_r_squared = adjusted_r_squared_values), aes(x = adjusted_r_squared)) +
  geom_histogram(binwidth = 0.05, fill = "steelblue", color = "black") +
  scale_x_continuous(limits = c(-1, 1), breaks = seq(-1, 1, by = 0.1)) +
  labs(title = "Distribution of Adjusted R-Squared for Significant Genes",
       x = "Adjusted R-Squared (-1 to 1)",
       y = "Frequency") +
  theme_minimal()
