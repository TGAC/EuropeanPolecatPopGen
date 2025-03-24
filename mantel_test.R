#this script carries out an isolation by distance on British polecat data using a mantel test. Input data are the geographic locations 
# (long/lat) of sample locations and genetic matrix which was calculated using the hamming distance in plink
#e.g. plink --vcf ${VCF} --allow-extra-chr --distance square --out  genetic_distance_UK_unrelated. 


setwd("/Users/shawr/Documents/PolecatWork/JAN2025/ibd/")
library(dplyr)
library(ggplot2)
library(geosphere)
library(ade4)

# Input: Data frame of individuals with lat/lon
coords <-read.csv("population_coords_unrelated.csv")


# Calculate pairwise geographic distances
geo_matrix_meters <- distm(coords[, c("Long", "Lat")], fun = distHaversine)
# Add row and column names (individual IDs)
rownames(geo_matrix_meters) <- coords$ID
colnames(geo_matrix_meters) <- coords$ID

# Convert the distances from meters to kilometers by dividing by 1000
geo_matrix_km <- geo_matrix_meters / 1000

# View the result (optional)
head(geo_matrix_km)

# If you'd like to save the result to a file
write.table(geo_matrix_km, "geo_distance_km.txt", sep = "\t", row.names = FALSE)

#load genetic matric from plink
genetic_matrix_pruned <- read.table("~/Documents/PolecatWork/JAN2025/ibd/genetic_distance_uk.dist_pruned.dist")
ids <- read.table("genetic_distance_uk.dist.id_pruned", header = FALSE)

# Add row and column names
rownames(genetic_matrix_pruned) <- ids$V1
colnames(genetic_matrix_pruned) <- ids$V1

# Perform Mantel test
mantel_result_pruned <- mantel.rtest(as.dist(genetic_matrix_pruned), as.dist(geo_matrix_km), nrepet = 999)

#plot correlation

# Flatten the upper triangle of both matrices into vectors
genetic_distances <- as.vector(as.dist(genetic_matrix_pruned))
geographic_distances <- as.vector(as.dist(geo_matrix_km))

# Extract individual IDs (from row/column names of the matrices)
ids <- colnames(genetic_matrix_pruned)  # Or rownames, depending on your matrix orientation

# Extract individual IDs from row/column names
ids <- rownames(genetic_matrix_pruned)  # Or colnames(genetic_matrix), depending on your matrix orientation

# Initialise empty lists to store pairwise values and IDs
genetic_distances <- c()
geographic_distances <- c()
Individual1 <- c()
Individual2 <- c()

# Loop through the matrix and collect pairwise comparisons
for (i in seq_len(nrow(genetic_matrix_pruned))) {
  for (j in seq_len(ncol(genetic_matrix_pruned))) {
    # Skip self-comparisons
    if (i != j) {
      genetic_distances <- c(genetic_distances, genetic_matrix_pruned[i, j])
      geographic_distances <- c(geographic_distances, geo_matrix_km[i, j])
      Individual1 <- c(Individual1, ids[i])
      Individual2 <- c(Individual2, ids[j])
    }
  }
}

# Combine into a data frame
data <- data.frame(
  GeneticDistance = genetic_distances,
  GeographicDistance = geographic_distances,
  Individual1 = Individual1,
  Individual2 = Individual2
)

# Step 2: Define the total number of SNPs in your dataset
total_snps <- 2941545 # Replace with the actual number of SNPs - taken from log file from run of plink

# Step 3: Normalize the Hamming distance (GeneticDistance)
data$NormalizedGeneticDistance <- data$GeneticDistance / total_snps

# Step 4: View the dataframe with the normalized Hamming distance
head(data)


# Create the plot
plot2 <- ggplot(data, aes(x = GeographicDistance, y = NormalizedGeneticDistance)) +
  geom_point(size = 4, alpha = 0.75, colour = "black", shape = 21, aes(fill = NormalizedGeneticDistance)) +
  
   #Add Individual1 labels (to the points)
  #geom_text(
   # aes(label = Individual2),  # Label points with Individual1 IDs
  #  size = 3, colour = "blue",  # Adjust size and colour of the labels
  #  nudge_y = 0.02, nudge_x = 0.02  # Nudge labels to avoid overlap
  #) +
  

  geom_smooth(method = "lm", colour = "black", alpha = 0.2) +  # Add trendline
  labs(
    x = "Geographic Distance (km)",
    y = "Genetic Distance (Hamming)",
    fill = "Genetic Distance"
  ) +
  theme(
    axis.text.x = element_text(face = "bold", colour = "black", size = 12),
    axis.text.y = element_text(face = "bold", size = 11, colour = "black"),
    axis.title = element_text(face = "bold", size = 14, colour = "black"),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black"),
    legend.position = "top",
    legend.text = element_text(size = 6, face = "bold"),
    legend.title = element_text(size = 11, face = "bold"),
    plot.background = element_blank(), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
  ) +
  scale_fill_continuous(
    low = "lightblue",
    high = "midnightblue",
    breaks = scales::breaks_pretty(n = 5)  # Automatically choose ~5 evenly spaced breaks
  )
