library(dplyr)

# e3/e3": This is the most common APOE genotype, considered the "neutral" or wild-type variant. It is associated with averattribute_4 risk for Alzheimer's disease and cardiovascular disorders.
# "e3/e4": This genotype carries one copy of the e3 allele and one copy of the e4 allele. It is associated with an increased risk of Alzheimer's disease, with individuals having approximately 3-4 times higher risk compared to those with e3/e3.
# "e4/e4": Homozygous for the e4 allele, this genotype is associated with the highest risk for Alzheimer's disease. Individuals with this genotype have about 12-15 times higher risk of developing Alzheimer's compared to those with e3/e3.
# "e2/e4": This genotype combines one protective allele (e2) with one risk allele (e4). Despite the presence of e2, individuals with this genotype still have an increased risk of Alzheimer's disease, similar to that of e3/e4 carriers.
# "e2/e3": This genotype is associated with a lower risk of Alzheimer's disease compared to e3/e3, due to the protective effect of the e2 allele34. However, it may be associated with an increased risk of type III hyperlipoproteinemia in some individuals

# Read the phenotype data
phenotype_data <- read.table("phenotype_file.txt", header = TRUE)

unique(phenotype_data$APOE)

phenotype_data <- phenotype_data %>%
  mutate(APOE_description = case_when(
    APOE == "e3/e3" ~ "Most common, averattribute_4 risk",
    APOE == "e3/e4" ~ "Increased risk for Alzheimer's",
    APOE == "e4/e4" ~ "Highest risk for Alzheimer's",
    APOE == "e2/e4" ~ "Mixed risk, still elevated",
    APOE == "e2/e3" ~ "Lower risk for Alzheimer's",
    TRUE ~ NA_character_
  ))

# Remove rows with NA in any column
df_clean <- phenotype_data %>% na.omit()

# Create a copy excluding ID and APOE columns
df_cluster <- df_clean %>% select(-ID, -APOE, -APOE4, -APOE_description)


# One-hot encoding for categorical variables
library(caret)

# Identify categorical columns
cat_cols <- c("STATUS", "SEX", "attribute_3")

# Perform one-hot encoding
df_encoded <- dummyVars(" ~ .", data = df_cluster) %>% predict(df_cluster)

# Convert back to data frame
df_encoded <- as.data.frame(df_encoded)

# Scale the numeric variables
df_scaled <- scale(df_encoded)



library(dplyr)
library(ggplot2)

# Function to calculate WCSS for a range of K values
calculate_wcss <- function(data, max_k) {
  sapply(1:max_k, function(k) {
    kmeans(data, centers = k, nstart = 10)$tot.withinss
  })
}

# Function to find the elbow point
find_elbow <- function(wcss) {
  n_points <- length(wcss)
  all_coords <- cbind(1:n_points, wcss)
  
  line_vec <- all_coords[n_points,] - all_coords[1,]
  line_vec_norm <- line_vec / sqrt(sum(line_vec^2))
  
  vec_from_first <- t(apply(all_coords, 1, function(coord) coord - all_coords[1,]))
  scalar_prod <- vec_from_first %*% line_vec_norm
  
  vec_from_line <- vec_from_first - scalar_prod %*% t(line_vec_norm)
  dist_from_line <- sqrt(rowSums(vec_from_line^2))
  
  which.max(dist_from_line)
}

# Calculate WCSS for a range of K values
max_k <- 40
wcss <- calculate_wcss(df_scaled, max_k)

# Find the optimal number of clusters
optimal_k <- find_elbow(wcss)

cat("The optimal number of clusters is:", optimal_k, "\n")

# Create the elbow plot with the optimal k indicated
elbow_data <- data.frame(k = 1:max_k, wcss = wcss)

ggplot(elbow_data, aes(x = k, y = wcss)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = optimal_k, linetype = "dashed", color = "red") +
  scale_x_continuous(breaks = 1:max_k) +
  labs(title = "Elbow Method for Optimal k",
       x = "Number of Clusters (k)",
       y = "Within-Cluster Sum of Squares") +
  theme_minimal() +
  annotate("text", x = optimal_k + 0.2, y = max(wcss), 
           label = paste("Optimal k =", optimal_k), 
           hjust = 0, vjust = 1)

# Perform clustering
library(cluster)
set.seed(123)


kmeans_result <- kmeans(df_scaled, centers = 5)

# Add cluster assignments back to the original dataframe
df_clean$Cluster <- kmeans_result$cluster

# Visualize the clustering results
library(ggplot2)

ggplot(df_clean, aes(x = attribute_4, y = attribute_2, color = factor(Cluster))) +
  geom_point() +
  facet_wrap(~APOE) +
  labs(title = "Clustering Results by APOE Genotype",
       x = "attribute_4", y = "attribute_2",
       color = "Cluster")


ggplot(df_clean, aes(x = attribute_4, y = attribute_3, color = factor(Cluster))) +
  geom_point() +
  facet_wrap(~APOE) +
  labs(title = "Clustering Results by APOE Genotype",
       x = "attribute_4", y = "attribute_3",
       color = "Cluster")



library(Rtsne)

set.seed(42)  # for reproducibility
tsne_result <- Rtsne(df_scaled, dims = 2, perplexity = 30, verbose = TRUE, max_iter = 500)

tsne_df <- data.frame(
  x = tsne_result$Y[,1],
  y = tsne_result$Y[,2],
  Cluster = factor(kmeans_result$cluster),
  APOE = df_clean$APOE
)

ggplot(tsne_df, aes(x = x, y = y, color = APOE, shape = Cluster)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(title = "t-SNE Visualization of Data",
       subtitle = "Colored by Cluster, Shaped by APOE Genotype",
       x = "t-SNE 1",
       y = "t-SNE 2") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "right")





ggplot(tsne_df, aes(x = x, y = y, color = Cluster)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(title = "t-SNE Visualization of Data",
       subtitle = "Colored by Cluster",
       x = "t-SNE 1",
       y = "t-SNE 2") +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "right")




ggplot(tsne_df, aes(x = x, y = y, color = APOE)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_minimal() +
  labs(title = "t-SNE Visualization of Data",
       subtitle = "Colored by APOE Genotype",
       x = "t-SNE 1",
       y = "t-SNE 2") +
  scale_color_brewer(palette = "Set1") +  # Changed to Set2 for different colors
  theme(legend.position = "right")




# Assuming your data frame is called 'df'
apoe_distribution <- table(phenotype_data$APOE)
apoe_percentattribute_4s <- prop.table(apoe_distribution) * 100
apoe_summary <- data.frame(
  Genotype = names(apoe_distribution),
  Count = as.vector(apoe_distribution),
  Percentattribute_4 = as.vector(apoe_percentattribute_4s)
)

# Print the summary table
print(apoe_summary)

ggplot(apoe_summary, aes(x = Genotype, y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = sprintf("%d\n(%.1f%%)", Count, Percentattribute_4)), 
            vjust = 0.5, size = 3) +
  labs(title = "Distribution of APOE Genotypes",
       x = "APOE Genotype",
       y = "Frequency") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))