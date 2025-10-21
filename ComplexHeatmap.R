library(Maaslin2)
require(Cairo)
library(ComplexHeatmap)
library(colorRamp2)
library(magick)
library(methods)
library(vegan)
library(dendsort)
library(openxlsx)
library(cluster)

#Define names for input files #Sample list has to be in the same order as metadata
#input_file <- paste("humann3_5_merged_pathabundance_unstratified_renorm.tsv")
input_file <- paste("SGB_species_inf_vs_uninf_PE.txt")
metadata_file <- paste("metadata_inf_vs_uninf_PE.txt")
#OTU_annotation <- paste("OTU_input_annotation.txt")
PWY_annot <- paste("PWY_input_annotation.txt")

#Import input table and metadata table
input_table <- read.table(input_file, sep="\t", header=T, row.names=1, stringsAsFactors=FALSE, comment.char="")
input_metadata <- read.table(metadata_file, sep="\t", header=T, row.names=1, stringsAsFactors=FALSE, comment.char="")
#input_OTU_annot <- read.table(OTU_annotation, sep="\t", header=T, row.names=1, stringsAsFactors=FALSE, comment.char="")
#input_KO_annot <- read.table(KO_annot, sep="\t", header=T, row.names=1, stringsAsFactors=FALSE, comment.char="")
input_PWY_annot <- read.table(PWY_annot, sep="\t", header=T, row.names=1, stringsAsFactors=FALSE, comment.char="")
#Removing to samples with "remove" flag in metadata file
specified_column <- c("remove")
specified_value <- c("remove")

######### Cluster by species
#### in the following lines, point towards the same files as above if species heatmap; or to a species table if doing a functional heatmap with species-based classification
#### In the metadata file, add a column name "remove", with "No" and "remove" values for each sample, all No's if you want to keep all

input_file_dendro_species <- paste("SGB_species_inf_vs_uninf_PE.txt")
metadata_file_dendro_species <- paste("metadata_inf_vs_uninf_PE.txt")
input_table_dendro_species <- read.table(input_file_dendro_species, sep="\t", header=T, row.names=1, stringsAsFactors=FALSE, comment.char="")
input_metadata_dendro_species <- read.table(metadata_file_dendro_species, sep="\t", header=T, row.names=1, stringsAsFactors=FALSE, comment.char="")
specified_column_dendro_species <- c("remove")
specified_value_dendro_species <- c("remove")
columns_to_remove_dendro_species <- names(input_table_dendro_species)[input_metadata_dendro_species[, specified_column_dendro_species] == specified_value_dendro_species]
rem_table_dendro_species <- input_table_dendro_species[, !names(input_table_dendro_species) %in% columns_to_remove_dendro_species]
clean_metadata_dendro_species <- subset(input_metadata_dendro_species, !remove %in% specified_value_dendro_species)
rounded_table_dendro_species <- round(rem_table_dendro_species)
col_sums_dendro_species <- colSums(rounded_table_dendro_species)
selected_cols_dendro_species <- col_sums_dendro_species > 0
nocol_filtered_dendro_species <- rounded_table_dendro_species[selected_cols_dendro_species]
clean_table_dendro_species <- as.data.frame(t(nocol_filtered_dendro_species))
dist_matrix_dendro_species <- vegdist(log(clean_table_dendro_species+1), method = "robust.aitchison")
hclust_obj_dendro_species <- hclust(dist_matrix_dendro_species, method = "ward.D2")
reord_hclust_obj_dendro_species <- dendsort(hclust_obj_dendro_species)
print(reord_hclust_obj_dendro_species)
plot(reord_hclust_obj_dendro_species)



# Identify columns to be removed based on the specified value in the metadata file
columns_to_remove <- names(input_table)[input_metadata[, specified_column] == specified_value]

# Filter out the specified columns from the sample file
rem_table <- input_table[, !names(input_table) %in% columns_to_remove]
clean_metadata <- subset(input_metadata, !remove %in% specified_value)

#Sum up rows to filter out empty KO/Taxonomy variables
rounded_table <- round(rem_table)

#Zero_col_remove
col_sums <- colSums(rounded_table)
selected_cols <- col_sums > 0
nocol_filtered <- rounded_table[selected_cols]

# #Zero_row_remove
# row_sums <- rowSums(nocol_filtered)
# input_filtered <- nocol_filtered[row_sums != 0, ]

#Average threshold filter to remove features in the final heatmap
row_means <- rowMeans(nocol_filtered)
input_filtered <- nocol_filtered[row_means > 10000, ]

# Convert to matrix and round up numbers
clean_table <- as.data.frame(t(nocol_filtered))

#Create dendogram object with external distance measures - this one is used for columns, copy-paste to create one for rows
dist_matrix <- vegdist(log(clean_table+1), method = "robust.aitchison")
hclust_obj <- hclust(dist_matrix, method = "ward.D2")
print(hclust_obj)
plot(hclust_obj)

reord_hclust_obj <- dendsort(hclust_obj)
print(reord_hclust_obj)
plot(reord_hclust_obj)

#Create dendogram object with external distance measures - this one is used for columns, copy-paste to create one for rows
rows_dist <- as.data.frame(t(clean_table))

# Count the number of positive values in each row
positive_counts <- rowSums(rows_dist > 0)

# Identify rows with two or fewer positive values
rows_to_remove <- positive_counts <= 10

# Subset the matrix, excluding the identified rows
filtered_matrix <- rows_dist[!rows_to_remove, ]

dist_matrix_rows <- vegdist(log(filtered_matrix+1), method = "jaccard")
dist_matrix_rows[is.na(dist_matrix_rows)] <- 0

# Step 5: Perform hierarchical clustering
row_clado <- hclust(as.dist(dist_matrix_rows), method = "ward.D2")

# Step 6: Plot dendrogram
plot(row_clado, main = "Row cluster", xlab = "Pathway_abundances")

reord_hclust_obj_rows <- dendsort(row_clado)
print(reord_hclust_obj_rows)
plot(reord_hclust_obj_rows)

# Reverse the order of the leaves
inverted_dend_rows <- reord_hclust_obj_rows
inverted_dend_rows$order <- rev(inverted_dend_rows$order)

# Plot the inverted dendrogram
plot(inverted_dend_rows)

nclust_rows <- 3

# Plot the dendrogram with highlighted clusters
rect.hclust(inverted_dend_rows, k = nclust_rows, border = "red")  # Highlight clusters with red rectangles

species_clusters <- cutree(inverted_dend_rows, nclust_rows)

# Extract labels and cluster belonging to data matrix
species_ids <- inverted_dend_rows$labels
species_cluster <- data.frame(sample_id = species_ids, cluster = species_clusters)

# Export as a tab-separated text file
write.table(species_cluster, file = "species_cluster.txt", sep = "\t", row.names = FALSE, col.names = FALSE)


#### WORKS BELOW

#Define the number of clusters to extract information
nclust <- 4
clusters <- cutree(reord_hclust_obj, nclust)

cluster_labels <- paste("Cluster", unique(clusters))
plot(reord_hclust_obj, labels = paste("Cluster", clusters))

# Extract labels and cluster belonging to data matrix
sample_ids <- reord_hclust_obj$labels
sample_cluster <- data.frame(sample_id = sample_ids, cluster = clusters)

#### Add the cluster annotation to the metadata file as an additional column

# Specify the column from the matrix to add to the metadata
column_to_add <- sample_cluster[rownames(clean_metadata), "cluster"]
clean_metadata$Functional_cluster <- column_to_add

# Specify the country variable
country_variable <- "Country"  # Replace with the actual variable name for country

# Specify the variables to summarize
variable1 <- "Infection_status"  # Replace with the actual variable name to summarize
variable2 <- "Functional_cluster"  # Replace with the actual variable name to summarize

# Stratify metadata by country variable
stratified_data <- split(clean_metadata, clean_metadata[[country_variable]])

# Create contingency tables for each country
contingency_tables <- lapply(stratified_data, function(subset) {
  table(subset[[variable1]], subset[[variable2]])
})

# Export contingency tables to individual text files
for (country in names(contingency_tables)) {
  table_name <- paste0("contingency_table_", country, ".txt")
  write.table(contingency_tables[[country]], file = table_name, quote = FALSE, sep = "\t")
}

# Highlight specified number of clusters
highlight_clusters <- 1:nclust
highlighted <- ifelse(clusters %in% highlight_clusters, clusters, NA)

# Export as a tab-separated text file
write.table(sample_cluster, file = "samples_cluster.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

# Plot the dendrogram with highlighted clusters
rect.hclust(reord_hclust_obj, k = nclust, border = "red")  # Highlight clusters with red rectangles

#Create dendogram object with external distance measures - this one is used for columns, copy-paste to create one for rows
#result <- spiec.easi(clean_table)

#diss_row_matrix
#bray_matrix <- vegdist(log(clean_table+1), method = "robust.aitchison")
#hclust_obj <- hclust(bray_matrix, method = "ward.D2")
#print(hclust_obj)
#plot(hclust_obj)

##reord_hclust_obj <- dendsort(hclust_obj)
#print(reord_hclust_obj)
#plot(reord_hclust_obj)

#Scale clean_table by column? to prepare for complex heatmap - uncomment below for col scaling
#scaled_clean_table <- scale(clean_table, scale = TRUE)
#scaled_clean_table <- as.data.frame(scaled_clean_table)

#Scale clean_table by row to prepare for complex heatmap - uncomment below for row scaling
clean_table <- t(clean_table)
scaled_clean_table <- scale(clean_table, scale = TRUE)
scaled_clean_table <- t(scaled_clean_table)
scaled_clean_table <- as.data.frame(scaled_clean_table)

#Create ComplexHeatmap compatible input
#heatmap_input <- as.matrix(t(scaled_clean_table))

#Create ComplexHeatmap compatible input
heatmap_input <- as.matrix(filtered_matrix)

#Create ComplexHeatmap compatible input using non-scaled table
#heatmap_input <- as.matrix(clean_table)

###not needed
#country_code <- c("Ivory_coast" = "darkorange", "Laos" = "blue4", "Pemba" = "black")

#Add country code to metadata file
#country_code <- c("Ivory_coast" = "1", "Laos" = "2", "Pemba" = "3")
#clean_metadata$color_country <- country_code[clean_metadata$Country]

#Extract country name
country_code_mapping <- clean_metadata$Country

#Extract infection status
infection_status_mapping <- clean_metadata$Infection_status

#EPG normalization (log10)
EPG <- clean_metadata$EPG
logEPG <- log10(EPG)
EPG_mapping <- logEPG
EPG_mapping[is.infinite(EPG_mapping)] <- 0

#OTU_annotation #first clean, and then extract
matrix_row_names <- rownames(filtered_matrix)
matrix_row_names <- as.character(matrix_row_names)

### Comment or uncomment one or the other

# # For OTU annotations
# metadata_column <- input_OTU_annot$OTU_ID # Replace "column_name" with the actual column name in your metadata
# metadata_column <- as.character(metadata_column)

#For PWY annotations
metadata_column <- input_PWY_annot$path_ID # Replace "column_name" with the actual column name in your metadata
metadata_column <- as.character(metadata_column)


matching_rows <- metadata_column %in% matrix_row_names

### Comment or uncomment one or the other

# # # For OTU annotations
# filtered_OTU_input <- input_OTU_annot[matching_rows, ]
# SGB_type_mapping <- filtered_OTU_input$Previously_described
# Kingdom_mapping <- filtered_OTU_input$Kingdom
# Phylum_mapping <- filtered_OTU_input$Phylum

# # For PWY annotations
filtered_PWY_input <- input_PWY_annot[matching_rows, ]
Superclass1 <- filtered_PWY_input$SUPERCLASS1
Superclass2 <- filtered_PWY_input$SUPERCLASS2

#Create annotation objects for ComplexHeatmap
#Country color
country_col <- c("Ivory_coast" = "#ff8c00", "Laos" = "blue4", "Pemba" = "black")

#Infection status
infect_col <- c("Infected" = "darkolivegreen4", "Not_infected"= "cornflowerblue")

#EPG counts
bottom_annotation = HeatmapAnnotation("log10(EPG)" = anno_points(EPG_mapping), annotation_height = unit(2, "cm"))

### Comment or uncomment one or the other

#SGB type col
SGB_col <- c("No" = "darkred", "Yes" = "darkolivegreen2", "Yes_species_level" = "white", "not_pres_SGB_DB" = "white")

#Create annotations
#Country color + #Infection status
#top_annotation = HeatmapAnnotation(Country = country_code_mapping, Infection_status = infection_status_mapping, col = list(Country = country_col, Infection_status = infect_col), simple_anno_size = unit(1, "cm"))

#NO Country color + #Infection status
top_annotation = HeatmapAnnotation(Infection_status = infection_status_mapping, col = list(Infection_status = infect_col), simple_anno_size = unit(1, "cm"))


### Comment or uncomment one or the other

# # OTU annotations
# right_annotation = rowAnnotation(SGB_type = SGB_type_mapping, Kingdom = Kingdom_mapping, Phylum = Phylum_mapping, col = list(SGB_type = SGB_col))

#PWY annotations
right_annotation = rowAnnotation(Superclass_1 = Superclass1, Superclass_2 = Superclass2)

# Customize the color scheme for the heatmap - works well for the functional data representation
col_fun <- colorRamp2(c(min(heatmap_input), quantile(heatmap_input, probs = 0.33), quantile(heatmap_input, probs = 0.66), max(heatmap_input)), colors = c("white", "darkgoldenrod1", "darkred", "black"))

# For a more granular color palette (e.g. taxonomy)

# Define the start and end colors
color_start <- "white"
color_middle1 <- "beige"
color_middle2 <- "darkgoldenrod1"
color_middle3 <- "darkorange"
color_middle4 <- "darkorange3"
color_end <- "darkred"

# Create the color gradient
#color_gradient <- colorRamp2(c(0, 10, 100, 1000, 10000, 100000), colors = c(color_start, color_middle1, color_middle2, color_middle3, color_middle4, color_end))
color_gradient <- colorRamp2(c(0, 1000, 10000, 100000), colors = c(color_start, color_middle2, color_middle4, color_end))

# Generate the gradient colors
num_colors <- 10000
gradient_colors <- color_gradient(seq(0, 450000, length.out = num_colors))

print(min(heatmap_input))
print(max(heatmap_input))

# # Define the color gradient range
# color_start <- "blue"
# color_middle <- "white"
# color_end <- "red"
# 
# # Get the minimum and maximum values in the matrix
# min_value <- min(heatmap_input)
# max_value <- max(heatmap_input)
# 
# # Define the levels for the gradient
# gradient_levels <- c(min(heatmap_input), mean(heatmap_input), max(heatmap_input))
# 
# # Create the color gradient function
# color_gradient <- colorRamp2(gradient_levels, colors = c(color_start, color_middle, color_end))
# 
# # Generate the gradient colors based on the matrix values
# gradient_colors <- color_gradient(heatmap_input)

# Create the initial Heatmap object
#heatmap_object <- Heatmap(heatmap_input, col = gradient_colors, name = "Species abundance", top_annotation = top_annotation, cluster_columns = reord_hclust_obj_dendro_species, cluster_rows = inverted_dend_rows, bottom_annotation = bottom_annotation, show_row_names = FALSE, right_annotation = right_annotation, show_column_names = FALSE, column_dend_height = unit(6, "cm"), row_dend_width = unit(6, "cm"))

# Create the initial Heatmap object
heatmap_object <- Heatmap(heatmap_input, col = gradient_colors, name = "Species abundance", top_annotation = top_annotation, cluster_columns = reord_hclust_obj_dendro_species, cluster_rows = inverted_dend_rows, bottom_annotation = bottom_annotation, show_row_names = FALSE, show_column_names = FALSE, column_dend_height = unit(6, "cm"), row_dend_width = unit(6, "cm"))

# Visualize the heatmap
draw(heatmap_object)

