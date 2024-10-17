library(tidyr)
library(ggplot2)
library(tidyverse)
library(mixOmics)
library(netOmics)
library(igraph)
library(gprofiler2)
library(paletteer)

# Set default image width
output_fold <- './output/plots/'

# Load data
load('./data/processed_data.RData')
cyto_scaled <- processed_data$cyto
prot_scaled <- processed_data$prot
rna_scaled <- processed_data$rna



#############################################################################
##################################   GRN   ##################################
#############################################################################
##### Shape data #####
cyto_wide <- cyto_scaled %>%
  dplyr::select(Participant, Condition, Cytokine, Value) %>%
  tidyr::pivot_wider(names_from = Cytokine, values_from = Value)
# Get rid of non numeric variables
cyto_matrix <- as.matrix(dplyr::select(cyto_wide, -Participant, -Condition))
condition_group <- cyto_wide$Condition

prot_wide <- prot_scaled %>%
  dplyr::select(Participant, Condition, Protein, Value) %>%
  tidyr::pivot_wider(names_from = Protein, values_from = Value)
# Get rid of non numeric variables
prot_matrix <- as.matrix(dplyr::select(prot_wide, -Participant, -Condition))

rna_wide <- rna_scaled %>%
  dplyr::select(Participant, Condition, RNA, Value) %>%
  tidyr::pivot_wider(names_from = RNA, values_from = Value)
# Get rid of non numeric variables
rna_matrix <- as.matrix(dplyr::select(rna_wide, -Participant, -Condition))


##### Build gene regulatory network #####
# grn only on genes
grn <- get_grn(rna_matrix)
vcount(grn)
ecount(grn)

png(paste0(output_fold, 'grn/grn_simple.png'),
    width = 3000, height = 3000, units = "px", pointsize = 72)
plot(grn, vertex.label = NA, vertex.color = "#a44200",
     vertex.size = 5,  edge.arrow.size = 0.5,
     main = "GRN on RNA data only")
dev.off()

non_connected_nodes <- V(grn)[degree(grn) == 0]
grn_clean <- delete_vertices(grn, non_connected_nodes)
png(paste0(output_fold, 'grn/grn_connected.png'),
    width = 3000, height = 3000, units = "px", pointsize = 72)
plot(grn_clean, vertex.label = NA, vertex.color = "#a44200",
     vertex.size = 5,  edge.arrow.size = 0.5,
     main = "GRN on RNA data after unconnected\nnodes deletion")
dev.off()

grn_deg <- degree(grn_clean)
png(paste0(output_fold, 'grn/grn_degrees.png'),
    width = 3000, height = 3000, units = "px", pointsize = 72)
hist(grn_deg,
     col = "#a44200", xlab = "Degree",
     ylab = "Frequency", main = "Degree Distribution of GRN")
dev.off()



#############################################################################
##################################   PPI   ##################################
#############################################################################
ppi <- readRDS("./data/supp/human_ppi.Rds")
protein_list <- unique(colnames(prot_matrix))
filtered_ppi <- ppi %>%
  filter(A %in% protein_list & B %in% protein_list)

ppi_prot <- graph_from_data_frame(filtered_ppi, directed = FALSE)
png(paste0(output_fold, 'grn/ppi.png'),
    width = 3000, height = 3000, units = "px", pointsize = 72)
plot(ppi_prot, vertex.label = NA, vertex.size = 5,
     edge.width = 2 , main = "PPI network", vertex.color = "#f2f1be")
dev.off()

ppi_deg <- degree(ppi_prot)
png(paste0(output_fold, 'grn/ppi_degrees.png'),
    width = 3000, height = 3000, units = "px", pointsize = 72)
hist(ppi_deg,
     col = "#f2f1be", xlab = "Degree",
     ylab = "Frequency", main = "Degree Distribution of PPI")
dev.off()



#############################################################################
#############################   Gene // prot   ##############################
#############################################################################
gene_to_prot <- readRDS("./data/supp/human_gene_to_coding_protein.Rds")
TF_to_gene <- readRDS("./data/supp/human_TF_to_targeted_gene.Rds")
gene_list <- V(grn_clean)$name


#### gene -> prot ####
gene_to_prot_c <- data.frame("from" = gene_to_prot$ENSEMBL, "to" = gene_to_prot$UNIPROT) %>%
  filter(from %in% gene_list & to %in% protein_list)

# Plot merged graph
merged_graph <- combine_layers(graph1 = grn_clean, graph2 = ppi_prot, interaction.df = gene_to_prot_c)
plot(merged_graph, vertex.label = NA, vertex.size = 5,
     edge.width = 2 , main = "Gene to prot connections")
# Colors
V(merged_graph)$color <- ifelse(V(merged_graph)$name %in% V(grn_clean)$name, "#a44200", "#f2f1be")
png(paste0(output_fold, 'grn/gene_to_prot.png'),
    width = 3000, height = 3000, units = "px", pointsize = 72)
plot(merged_graph, vertex.label = NA, vertex.size = 5,
     edge.width = 2 , main = "Gene to prot connections")
dev.off()

# Number of prot/gene connections
as.data.frame(as_data_frame(merged_graph, what = "edges")) %>%
  filter((from %in% gene_list & to %in% protein_list) |
           (to %in% gene_list & from %in% protein_list))


#### TF -> gene ####
TF_to_gene_c <- data.frame("from" = TF_to_gene$UNIPROT, "to" = TF_to_gene$ENSEMBL) %>%
  filter(from %in% protein_list & to %in% gene_list)

# Plot merged graph
merged_graph <- combine_layers(graph1 = grn_clean, graph2 = ppi_prot, interaction.df = TF_to_gene_c)
plot(merged_graph, vertex.label = NA, vertex.size = 5,
     edge.width = 2 , main = "TF to gene connections")
# Colors
V(merged_graph)$color <- ifelse(V(merged_graph)$name %in% V(grn_clean)$name, "#a44200", "#f2f1be")
png(paste0(output_fold, 'grn/TF_to_gene.png'),
    width = 3000, height = 3000, units = "px", pointsize = 72)
plot(merged_graph, vertex.label = NA, vertex.size = 5,
     edge.width = 2 , main = "TF to gene connections")
dev.off()

# Number of prot/gene connections
as.data.frame(as_data_frame(merged_graph, what = "edges")) %>%
  filter((from %in% gene_list & to %in% protein_list) |
           (to %in% gene_list & from %in% protein_list))


##### Global interactions #####
gene_prot_connect <- rbind(
  data.frame("from" = gene_to_prot$UNIPROT, "to" = gene_to_prot$ENSEMBL),
  data.frame("from" = TF_to_gene$UNIPROT, "to" = TF_to_gene$ENSEMBL)) %>%
  filter(from %in% protein_list & to %in% gene_list)

# Plot merged graph
merged_graph <- combine_layers(graph1 = grn_clean, graph2 = ppi_prot, interaction.df = gene_prot_connect)
plot(merged_graph, vertex.label = NA, vertex.size = 5,
     edge.width = 2 , main = "Connected gene and protein networks")
# Colors
V(merged_graph)$color <- ifelse(V(merged_graph)$name %in% V(grn_clean)$name, "#a44200", "#f2f1be")
png(paste0(output_fold, 'grn/grn_ppi.png'),
    width = 3000, height = 3000, units = "px", pointsize = 72)
plot(merged_graph, vertex.label = NA, vertex.size = 5,
     edge.width = 2 , main = "Connected gene and protein networks")
dev.off()


# Analyse de modularité
save(merged_graph, file = "./data/merged_graph.RData")
clusters_inter <- cluster_optimal(merged_graph)
plot(clusters_inter, merged_graph)


##### Random walk #####
random_walk <- random_walk(merged_graph, start = "P17676", steps = 500)
visited_nodes <- unique(random_walk)$name
length(visited_nodes)
V(merged_graph)$color <- ifelse(V(merged_graph)$name %in% visited_nodes, "#69150f",  "darkgrey")
png(paste0(output_fold, 'grn/random_walk.png'),
    width = 3000, height = 3000, units = "px", pointsize = 72)
plot(merged_graph, vertex.label = NA, vertex.size = 5,
     edge.width = 2, main = "Random Walk from Node P17676\n(500 steps)")
dev.off()


##### Enrichment #####
# Select main network
components <- clusters(merged_graph)
membership <- components$membership
# Find out how many nodes are in each component
component_sizes <- components$csize
# Check the number of components and their sizes
print(component_sizes)
# Extract the biggest
largest_component_id <- which.max(component_sizes)
largest_subgraph <- induced_subgraph(merged_graph, which(membership == largest_component_id))

# Node names
net_list <- V(largest_subgraph)$name
gostres <- gost(query = net_list, organism = "hsapiens")

# Analysis
png(paste0(output_fold, 'grn/enrichment.png'),
    width = 3000, height = 3000, units = "px", pointsize = 1000)
gostplot(gostres, capped = FALSE, interactive = FALSE,
  pal = paletteer_c("ggthemes::Orange-Gold", 11))
dev.off()

go_cc_results <- gostres$result[gostres$result$source == "GO:CC", ] %>%
  arrange(intersection_size)
go_bp_results <- gostres$result[gostres$result$source == "GO:BP", ] %>%
  arrange(intersection_size)
go_bp_results$short_term_name <- substr(go_bp_results$term_name, 1, 50)

ggplot(go_cc_results[1:30,], aes(x = reorder(term_name, -intersection_size), y = intersection_size)) +
  geom_bar(stat = "identity", fill = "#f2f1be") +  # Change color as needed
  coord_flip() +  # Flip coordinates for better readability
  xlab("GO Term (Cellular Component)") +
  ylab("Intersection Size (Number of Genes)") +
  ggtitle("Top 30 Gene Overlaps\nwith GO:CC Terms") +
  theme_minimal(base_size = 12) +
  ylim(c(0, 20))

ggplot(go_bp_results[1:30,], aes(x = reorder(short_term_name, -intersection_size), y = intersection_size)) +
  geom_bar(stat = "identity", fill = "#d58937") +  # Change color as needed
  coord_flip() +  # Flip coordinates for better readability
  xlab("GO Term (Biological Process)") +
  ylab("Intersection Size (Number of Genes)") +
  ggtitle("Top 30 Gene Overlaps\nwith GO:BP Terms") +
  theme_minimal(base_size = 12) +
  ylim(c(0, 20))








