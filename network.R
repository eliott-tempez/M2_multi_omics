library(tidyr)
library(ggplot2)
library(tidyverse)
library(mixOmics)
library(netOmics)
library(igraph)
library(timeOmics)

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
     edge.width = 2 , main = "PPI network", vertex.color = "#a44200")
dev.off()

ppi_deg <- degree(ppi_prot)
png(paste0(output_fold, 'grn/ppi_degrees.png'),
    width = 3000, height = 3000, units = "px", pointsize = 72)
hist(ppi_deg,
     col = "#a44200", xlab = "Degree",
     ylab = "Frequency", main = "Degree Distribution of PPI")
dev.off()



#############################################################################
#############################   Gene // prot   ##############################
#############################################################################
gene_to_prot <- readRDS("./data/supp/human_gene_to_coding_protein.Rds")
TF_to_gene <- readRDS("./data/supp/human_TF_to_targeted_gene.Rds")











