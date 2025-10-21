# Load necessary libraries
library(SpiecEasi)
library(igraph)
library(ggplot2)
library(tools)
library(RCy3)

#Output folder
output_folder <- '20251021_output_comparisons_200_species'
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

# Paths to the input files
abundance_file_path_1 <- 'inputs/IC_infected.txt'
abundance_file_path_2 <- 'inputs/IC_uninfected.txt'
abundance_file_path_3 <- 'inputs/LA_infected.txt'
abundance_file_path_4 <- 'inputs/LA_uninfected.txt'
abundance_file_path_5 <- 'inputs/PE_infected.txt'
abundance_file_path_6 <- 'inputs/PE_uninfected.txt'
taxa_file_path <- 'highlight_taxa_list.txt'  # Path to the file with taxa to highlight

# Extract base names from file paths for naming output files
abundance_file_name_1 <- file_path_sans_ext(basename(abundance_file_path_1))
abundance_file_name_2 <- file_path_sans_ext(basename(abundance_file_path_2))
abundance_file_name_3 <- file_path_sans_ext(basename(abundance_file_path_3))
abundance_file_name_4 <- file_path_sans_ext(basename(abundance_file_path_4))
abundance_file_name_5 <- file_path_sans_ext(basename(abundance_file_path_5))
abundance_file_name_6 <- file_path_sans_ext(basename(abundance_file_path_6))

# Read the abundance tables
abundance_table_1 <- read.table(abundance_file_path_1, header = TRUE, sep = "\t", row.names = 1)
abundance_table_1 <- t(abundance_table_1)

abundance_table_2 <- read.table(abundance_file_path_2, header = TRUE, sep = "\t", row.names = 1)
abundance_table_2 <- t(abundance_table_2)

abundance_table_3 <- read.table(abundance_file_path_3, header = TRUE, sep = "\t", row.names = 1)
abundance_table_3 <- t(abundance_table_3)

abundance_table_4 <- read.table(abundance_file_path_4, header = TRUE, sep = "\t", row.names = 1)
abundance_table_4 <- t(abundance_table_4)

abundance_table_5 <- read.table(abundance_file_path_5, header = TRUE, sep = "\t", row.names = 1)
abundance_table_5 <- t(abundance_table_5)

abundance_table_6 <- read.table(abundance_file_path_6, header = TRUE, sep = "\t", row.names = 1)
abundance_table_6 <- t(abundance_table_6)

# Keep only the top 50 species with the highest average abundance
avg_abundance_1 <- colMeans(abundance_table_1)
top_50_species_1 <- names(sort(avg_abundance_1, decreasing = TRUE))[1:200]
abundance_table_1 <- abundance_table_1[, top_50_species_1, drop = FALSE]

avg_abundance_2 <- colMeans(abundance_table_2)
top_50_species_2 <- names(sort(avg_abundance_2, decreasing = TRUE))[1:200]
abundance_table_2 <- abundance_table_2[, top_50_species_2, drop = FALSE]

avg_abundance_3 <- colMeans(abundance_table_3)
top_50_species_3 <- names(sort(avg_abundance_3, decreasing = TRUE))[1:200]
abundance_table_3 <- abundance_table_3[, top_50_species_3, drop = FALSE]

avg_abundance_4 <- colMeans(abundance_table_4)
top_50_species_4 <- names(sort(avg_abundance_4, decreasing = TRUE))[1:200]
abundance_table_4 <- abundance_table_4[, top_50_species_4, drop = FALSE]

avg_abundance_5 <- colMeans(abundance_table_5)
top_50_species_5 <- names(sort(avg_abundance_5, decreasing = TRUE))[1:200]
abundance_table_5 <- abundance_table_5[, top_50_species_5, drop = FALSE]

avg_abundance_6 <- colMeans(abundance_table_6)
top_50_species_6 <- names(sort(avg_abundance_6, decreasing = TRUE))[1:200]
abundance_table_6 <- abundance_table_6[, top_50_species_6, drop = FALSE]

# Read the list of taxa to highlight
highlight_taxa_list <- readLines(taxa_file_path)
highlight_taxa_list <- highlight_taxa_list[highlight_taxa_list != ""]  # Remove any empty lines

# Run SpiecEasi
se_1 <- spiec.easi(as.matrix(abundance_table_1), method = 'mb', nlambda = 20, lambda.min.ratio = 1e-2, pulsar.params = list(rep.num = 100))
se_2 <- spiec.easi(as.matrix(abundance_table_2), method = 'mb', nlambda = 20, lambda.min.ratio = 1e-2, pulsar.params = list(rep.num = 100))
se_3 <- spiec.easi(as.matrix(abundance_table_3), method = 'mb', nlambda = 20, lambda.min.ratio = 1e-2, pulsar.params = list(rep.num = 100))
se_4 <- spiec.easi(as.matrix(abundance_table_4), method = 'mb', nlambda = 20, lambda.min.ratio = 1e-2, pulsar.params = list(rep.num = 100))
se_5 <- spiec.easi(as.matrix(abundance_table_5), method = 'mb', nlambda = 20, lambda.min.ratio = 1e-2, pulsar.params = list(rep.num = 100))
se_6 <- spiec.easi(as.matrix(abundance_table_6), method = 'mb', nlambda = 20, lambda.min.ratio = 1e-2, pulsar.params = list(rep.num = 100))

# Generate networks
adj_1 <- getRefit(se_1)
rownames(adj_1) <- colnames(adj_1) <- colnames(abundance_table_1)
g1 <- graph.adjacency(adj_1, mode = "undirected", diag = FALSE)

adj_2 <- getRefit(se_2)
rownames(adj_2) <- colnames(adj_2) <- colnames(abundance_table_2)
g2 <- graph.adjacency(adj_2, mode = "undirected", diag = FALSE)

adj_3 <- getRefit(se_3)
rownames(adj_3) <- colnames(adj_3) <- colnames(abundance_table_3)
g3 <- graph.adjacency(adj_3, mode = "undirected", diag = FALSE)

adj_4 <- getRefit(se_4)
rownames(adj_4) <- colnames(adj_4) <- colnames(abundance_table_4)
g4 <- graph.adjacency(adj_4, mode = "undirected", diag = FALSE)

adj_5 <- getRefit(se_5)
rownames(adj_5) <- colnames(adj_5) <- colnames(abundance_table_5)
g5 <- graph.adjacency(adj_5, mode = "undirected", diag = FALSE)

adj_6 <- getRefit(se_6)
rownames(adj_6) <- colnames(adj_6) <- colnames(abundance_table_6)
g6 <- graph.adjacency(adj_6, mode = "undirected", diag = FALSE)

# Combine the vertex names from both graphs to ensure consistent numbering
all_vertices <- unique(c(V(g1)$name, V(g2)$name, V(g3)$name, V(g4)$name, V(g5)$name, V(g6)$name))
vertex_legend <- data.frame(Number = 1:length(all_vertices), Label = all_vertices, stringsAsFactors = FALSE)

# Create a mapping from vertex names to numbers
vertex_mapping <- setNames(vertex_legend$Number, vertex_legend$Label)

# Apply the consistent numbering to both graphs
V(g1)$number <- vertex_mapping[V(g1)$name]
V(g2)$number <- vertex_mapping[V(g2)$name]
V(g3)$number <- vertex_mapping[V(g3)$name]
V(g4)$number <- vertex_mapping[V(g4)$name]
V(g5)$number <- vertex_mapping[V(g5)$name]
V(g6)$number <- vertex_mapping[V(g6)$name]

# Replace underscores with spaces in vertex names
V(g1)$name <- gsub("_", " ", V(g1)$name)
V(g2)$name <- gsub("_", " ", V(g2)$name)
V(g3)$name <- gsub("_", " ", V(g3)$name)
V(g4)$name <- gsub("_", " ", V(g4)$name)
V(g5)$name <- gsub("_", " ", V(g5)$name)
V(g6)$name <- gsub("_", " ", V(g6)$name)

# Highlight specified taxa in individual networks
V(g1)$color <- ifelse(V(g1)$name %in% gsub("_", " ", highlight_taxa_list), 'red', 'blue')
V(g2)$color <- ifelse(V(g2)$name %in% gsub("_", " ", highlight_taxa_list), 'red', 'blue')
V(g3)$color <- ifelse(V(g3)$name %in% gsub("_", " ", highlight_taxa_list), 'red', 'blue')
V(g4)$color <- ifelse(V(g4)$name %in% gsub("_", " ", highlight_taxa_list), 'red', 'blue')
V(g5)$color <- ifelse(V(g5)$name %in% gsub("_", " ", highlight_taxa_list), 'red', 'blue')
V(g6)$color <- ifelse(V(g6)$name %in% gsub("_", " ", highlight_taxa_list), 'red', 'blue')

# Combine the 6 networks
combined_edges <- unique(rbind(as.data.frame(get.edgelist(g1)), as.data.frame(get.edgelist(g2)), as.data.frame(get.edgelist(g3)), as.data.frame(get.edgelist(g4)), as.data.frame(get.edgelist(g5)), as.data.frame(get.edgelist(g6))))
combined_graph <- graph_from_data_frame(combined_edges, directed = FALSE)

# Replace underscores with spaces in combined network vertex names
V(combined_graph)$name <- gsub("_", " ", V(combined_graph)$name)

# Apply the consistent numbering to the combined graph
V(combined_graph)$number <- vertex_mapping[V(combined_graph)$name]

# Highlight specified taxa in combined network
V(combined_graph)$color <- ifelse(V(combined_graph)$name %in% gsub("_", " ", highlight_taxa_list), 'red', 'blue')

# Highlight edges present in both, only in g1, and only in g2
#E(combined_graph)$color <- "gray"
#E(combined_graph)[which(E(combined_graph) %in% E(g1) & E(combined_graph) %in% E(g2))]$color <- "green" # Present in both
#E(combined_graph)[which(E(combined_graph) %in% E(g1) & !(E(combined_graph) %in% E(g2)))]$color <- "blue"  # Only in g1
#E(combined_graph)[which(!(E(combined_graph) %in% E(g1)) & E(combined_graph) %in% E(g2))]$color <- "orange"  # Only in g2

# Calculate network properties
properties_g1 <- list(
  num_nodes = vcount(g1),
  num_edges = ecount(g1),
  density = graph.density(g1),
  average_path_length = average.path.length(g1, directed = FALSE),
  diameter = diameter(g1, directed = FALSE),
  clustering_coefficient = transitivity(g1, type = "global")
)

properties_g2 <- list(
  num_nodes = vcount(g2),
  num_edges = ecount(g2),
  density = graph.density(g2),
  average_path_length = average.path.length(g2, directed = FALSE),
  diameter = diameter(g2, directed = FALSE),
  clustering_coefficient = transitivity(g2, type = "global")
)

properties_g3 <- list(
  num_nodes = vcount(g3),
  num_edges = ecount(g3),
  density = graph.density(g3),
  average_path_length = average.path.length(g3, directed = FALSE),
  diameter = diameter(g3, directed = FALSE),
  clustering_coefficient = transitivity(g3, type = "global")
)

properties_g4 <- list(
  num_nodes = vcount(g4),
  num_edges = ecount(g4),
  density = graph.density(g4),
  average_path_length = average.path.length(g4, directed = FALSE),
  diameter = diameter(g4, directed = FALSE),
  clustering_coefficient = transitivity(g4, type = "global")
)

properties_g5 <- list(
  num_nodes = vcount(g5),
  num_edges = ecount(g5),
  density = graph.density(g5),
  average_path_length = average.path.length(g5, directed = FALSE),
  diameter = diameter(g5, directed = FALSE),
  clustering_coefficient = transitivity(g5, type = "global")
)

properties_g6 <- list(
  num_nodes = vcount(g6),
  num_edges = ecount(g6),
  density = graph.density(g6),
  average_path_length = average.path.length(g6, directed = FALSE),
  diameter = diameter(g6, directed = FALSE),
  clustering_coefficient = transitivity(g6, type = "global")
)

# Compare centrality measures
centrality_g1 <- data.frame(
  species = V(g1)$name,
  degree = degree(g1),
  betweenness = betweenness(g1),
  closeness = closeness(g1),
  eigenvector = eigen_centrality(g1)$vector
)

centrality_g2 <- data.frame(
  species = V(g2)$name,
  degree = degree(g2),
  betweenness = betweenness(g2),
  closeness = closeness(g2),
  eigenvector = eigen_centrality(g2)$vector
)

centrality_g3 <- data.frame(
  species = V(g3)$name,
  degree = degree(g3),
  betweenness = betweenness(g3),
  closeness = closeness(g3),
  eigenvector = eigen_centrality(g3)$vector
)

centrality_g4 <- data.frame(
  species = V(g4)$name,
  degree = degree(g4),
  betweenness = betweenness(g4),
  closeness = closeness(g4),
  eigenvector = eigen_centrality(g4)$vector
)

centrality_g5 <- data.frame(
  species = V(g5)$name,
  degree = degree(g5),
  betweenness = betweenness(g5),
  closeness = closeness(g5),
  eigenvector = eigen_centrality(g5)$vector
)

centrality_g6 <- data.frame(
  species = V(g6)$name,
  degree = degree(g6),
  betweenness = betweenness(g6),
  closeness = closeness(g6),
  eigenvector = eigen_centrality(g6)$vector
)

#centrality_comparison <- merge(centrality_g1, centrality_g2, centrality_g3, centrality_g4, centrality_g5, centrality_g6, by = "species", suffixes = c("_g1", "_g2", "_g3", "_g4", "_g5", "_g6"))

# Compare edges
edges_g1 <- as.data.frame(get.edgelist(g1))
edges_g2 <- as.data.frame(get.edgelist(g2))
edges_g3 <- as.data.frame(get.edgelist(g3))
edges_g4 <- as.data.frame(get.edgelist(g4))
edges_g5 <- as.data.frame(get.edgelist(g5))
edges_g6 <- as.data.frame(get.edgelist(g6))
#common_edges <- merge(edges_g1, edges_g2, edges_g3, edges_g4, edges_g5, edges_g6, by = c("V1", "V2"))
#num_common_edges <- nrow(common_edges)

# Print results
print("Network Properties for Network 1:")
print(properties_g1)

print("Network Properties for Network 2:")
print(properties_g2)

print("Network Properties for Network 3:")
print(properties_g3)

print("Network Properties for Network 4:")
print(properties_g4)

print("Network Properties for Network 5:")
print(properties_g5)

print("Network Properties for Network 6:")
print(properties_g6)

print("Centrality Measures Comparison:")
#print(centrality_comparison)

#print(paste("Number of common edges:", num_common_edges))

# Save the centrality measures and network plots
#write.csv(centrality_comparison, file.path(output_folder, paste0('centrality_comparison_', abundance_file_name_1, '_', abundance_file_name_2, '_', abundance_file_name_3, '_', abundance_file_name_4, '_', abundance_file_name_5, '_', abundance_file_name_6, '.csv')), row.names = FALSE)

# Function to plot network with node numbers and corresponding legend
plot_network <- function(graph, output_path_tiff, output_path_svg) {
  coords <- layout_with_fr(graph)  # Fruchterman-Reingold layout
  
  # Plot to TIFF
  tiff(output_path_tiff, width = 1200, height = 1200)
  plot(graph, vertex.label = V(graph)$number, vertex.size = 13, vertex.label.color="white", vertex.label.cex = 4, edge.width = 2,
       vertex.frame.color = NA, edge.color = E(graph)$color, layout = coords)
  dev.off()
  
  # Plot to SVG
  svg(output_path_svg, width = 10, height = 10)
  plot(graph, vertex.label = V(graph)$number, vertex.size = 13, vertex.label.color="white", vertex.label.cex = 4, edge.width = 2,
       vertex.frame.color = NA, edge.color = E(graph)$color, layout = coords)
  dev.off()
}

# Plot and save individual and combined networks
plot_network(g1, 
             file.path(output_folder, paste0('network_plot_', abundance_file_name_1, '.tiff')),
             file.path(output_folder, paste0('network_plot_', abundance_file_name_1, '.svg')))

plot_network(g2, 
             file.path(output_folder, paste0('network_plot_', abundance_file_name_2, '.tiff')),
             file.path(output_folder, paste0('network_plot_', abundance_file_name_2, '.svg')))

plot_network(g3, 
             file.path(output_folder, paste0('network_plot_', abundance_file_name_3, '.tiff')),
             file.path(output_folder, paste0('network_plot_', abundance_file_name_3, '.svg')))

plot_network(g4, 
             file.path(output_folder, paste0('network_plot_', abundance_file_name_4, '.tiff')),
             file.path(output_folder, paste0('network_plot_', abundance_file_name_4, '.svg')))

plot_network(g5, 
             file.path(output_folder, paste0('network_plot_', abundance_file_name_5, '.tiff')),
             file.path(output_folder, paste0('network_plot_', abundance_file_name_5, '.svg')))

plot_network(g6, 
             file.path(output_folder, paste0('network_plot_', abundance_file_name_6, '.tiff')),
             file.path(output_folder, paste0('network_plot_', abundance_file_name_6, '.svg')))

plot_network(combined_graph, 
             file.path(output_folder, paste0('combined_network_plot_', abundance_file_name_1, '_', abundance_file_name_2, '_', abundance_file_name_3, '_', abundance_file_name_4, '_', abundance_file_name_5, '_', abundance_file_name_6, '.tiff')),
             file.path(output_folder, paste0('combined_network_plot_', abundance_file_name_1, '_', abundance_file_name_2, '_', abundance_file_name_3, '_', abundance_file_name_4, '_', abundance_file_name_5, '_', abundance_file_name_6, '.svg')))

# Save the unified legend
write.csv(vertex_legend, file.path(output_folder, 'unified_legend.csv'), row.names = FALSE)

# Save network properties
properties_g1_df <- as.data.frame(t(unlist(properties_g1)))
properties_g2_df <- as.data.frame(t(unlist(properties_g2)))
properties_g3_df <- as.data.frame(t(unlist(properties_g3)))
properties_g4_df <- as.data.frame(t(unlist(properties_g4)))
properties_g5_df <- as.data.frame(t(unlist(properties_g5)))
properties_g6_df <- as.data.frame(t(unlist(properties_g6)))
write.csv(properties_g1_df, file.path(output_folder, paste0('network_properties_', abundance_file_name_1, '.csv')), row.names = FALSE)
write.csv(properties_g2_df, file.path(output_folder, paste0('network_properties_', abundance_file_name_2, '.csv')), row.names = FALSE)
write.csv(properties_g3_df, file.path(output_folder, paste0('network_properties_', abundance_file_name_3, '.csv')), row.names = FALSE)
write.csv(properties_g4_df, file.path(output_folder, paste0('network_properties_', abundance_file_name_4, '.csv')), row.names = FALSE)
write.csv(properties_g5_df, file.path(output_folder, paste0('network_properties_', abundance_file_name_5, '.csv')), row.names = FALSE)
write.csv(properties_g6_df, file.path(output_folder, paste0('network_properties_', abundance_file_name_6, '.csv')), row.names = FALSE)

# Function to export network to Cytoscape
# export_to_cytoscape <- function(graph, network_name) {
#   # Ensure Cytoscape is running and accessible
#   cytoscapePing()
#   
#   # Create a data frame for nodes
#   nodes <- data.frame(id = V(graph)$name, label = V(graph)$name, color = V(graph)$color, number = V(graph)$number, stringsAsFactors = FALSE)
#   
#   # Create a data frame for edges
#   edges <- as.data.frame(get.edgelist(graph), stringsAsFactors = FALSE)
#   colnames(edges) <- c("source", "target")
#   edges$interaction <- "pp"
#   edges$color <- E(graph)$color
#   
#   # Create the network in Cytoscape
#   createNetworkFromDataFrames(nodes, edges, title = network_name, collection = "Microbiome Networks")
# }
# 
# # Export individual and combined networks to Cytoscape
# export_to_cytoscape(g1, paste0('Network_', abundance_file_name_1))
# export_to_cytoscape(g2, paste0('Network_', abundance_file_name_2))
# export_to_cytoscape(combined_graph, paste0('Combined_Network_', abundance_file_name_1, '_', abundance_file_name_2))
