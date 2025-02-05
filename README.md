# APOE Phenotype Analysis and Clustering

## Overview
This project analyzes the distribution of APOE genotypes and their association with Alzheimer's disease risk. It performs data preprocessing, one-hot encoding, scaling, clustering, and visualization using various techniques including k-means and t-SNE.

## Dependencies
The following R packages are required:

```r
install.packages(c("dplyr", "ggplot2", "caret", "cluster", "Rtsne"))
```

## Data Processing Steps
1. **Read Phenotype Data**  
   - Loads the dataset containing APOE genotypes and other attributes.
2. **Genotype Risk Classification**  
   - Assigns risk categories to different APOE genotypes.
3. **Data Cleaning & Encoding**  
   - Removes rows with missing values.
   - Performs one-hot encoding on categorical variables.
   - Scales numeric variables.

## Clustering Analysis
1. **Optimal Cluster Determination**  
   - Computes Within-Cluster Sum of Squares (WCSS).
   - Identifies the optimal number of clusters using the elbow method.
2. **K-Means Clustering**  
   - Assigns clusters to each sample based on scaled attributes.
   - Visualizes clustering results using scatter plots.

## t-SNE Visualization
- Applies t-SNE for dimensionality reduction.
- Generates visualizations based on clusters and APOE genotypes.

## Genotype Distribution
- Computes frequency and percentage of different APOE genotypes.
- Generates a bar chart for genotype distribution.

## Running the Code
Ensure the dataset `phenotype_file.txt` is available in the working directory, then execute the script in R.

```r
Rscript script_name.R
```

## Output
- Clustering assignments.
- Elbow plot for optimal cluster selection.
- t-SNE visualizations.
- Bar chart for APOE genotype distribution.

## Contact
For questions or contributions, please reach out via GitHub or email.
