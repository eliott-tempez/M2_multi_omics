library(tidyr)
library(ggplot2)
library(tidyverse)


#############################################################################
########################   Présentation du dataset   ########################
#############################################################################
data <- readRDS('./data/data_pregnancy.Rds')
# Number of tables in data
ls(data)

# Shape of each table
dim(data$RNA)
dim(data$protein)
dim(data$cyto)
dim(data$sample_info)
length(unique(data$sample_info$Y))

# On a 68 samples : 17 participantes et 4 conditions (3 trimestres + post-partum)
# Pour 3 types de données : ARN, protéines, cytokines
# Respectivement : 23855, 2329, 62 features


#############################################################################
#############################   Préprocessing   #############################
#############################################################################

##### CYTOKINES #####
# Transform data
cyto <- data$cyto %>%
  # Add a new column splitting the rownames into "Participant" and "Condition"
  mutate(Participant = sub("_.*", "", rownames(data$cyto)),
         Condition = sub(".*_", "", rownames(data$cyto))) %>%
  # Reshape from wide to long format
  pivot_longer(cols = -c(Participant, Condition), # Pivot all columns except Participant and Condition
               names_to = "Cytokine",             # The names of the columns will go into "Cytokine"
               values_to = "Value") %>%           # The corresponding values will go into "Value"
  arrange(Condition)

ggplot(cyto, aes(x = Cytokine, y = Value, col = Condition)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Cytokine Levels", 
       x = "Cytokine", 
       y = "Expression level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


ggplot(cyto, aes(x = Condition, y = Value, col = Condition)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Cytokine Levels by Condition", 
       x = "Time of prelevement", 
       y = "Expression level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggplot(cyto, aes(x = Value, fill = Condition)) +
  geom_histogram() +
  facet_wrap(vars(Condition)) +
  labs(x = "Cytokine expression") +
  theme(legend.position = "none")


##### PROTEINS #####
# Transform data
prot <- data$protein %>%
  # Add a new column splitting the rownames into "Participant" and "Condition"
  mutate(Participant = sub("_.*", "", rownames(data$cyto)),
         Condition = sub(".*_", "", rownames(data$cyto))) %>%
  # Reshape from wide to long format
  pivot_longer(cols = -c(Participant, Condition), # Pivot all columns except Participant and Condition
               names_to = "Protein",             # The names of the columns will go into "Cytokine"
               values_to = "Value") %>%           # The corresponding values will go into "Value"
  arrange(Condition)


ggplot(prot, aes(x = Condition, y = Value, col = Condition)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Protein Levels by Condition", 
       x = "Time of prelevement", 
       y = "Expression level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggplot(prot, aes(x = Value, fill = Condition)) +
  geom_histogram() +
  facet_wrap(vars(Condition)) +
  labs(x = "Protein expression") +
  theme(legend.position = "none")


##### RNA #####
# Transform data
rna <- data$RNA %>%
  # Add a new column splitting the rownames into "Participant" and "Condition"
  mutate(Participant = sub("_.*", "", rownames(data$cyto)),
         Condition = sub(".*_", "", rownames(data$cyto))) %>%
  # Reshape from wide to long format
  pivot_longer(cols = -c(Participant, Condition), # Pivot all columns except Participant and Condition
               names_to = "RNA",             # The names of the columns will go into "Cytokine"
               values_to = "Value") %>%           # The corresponding values will go into "Value"
  arrange(Condition)


ggplot(rna, aes(x = Condition, y = Value, col = Condition)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "RNA Levels by Condition", 
       x = "Time of prelevement", 
       y = "Expression level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggplot(rna, aes(x = Value, fill = Condition)) +
  geom_histogram() +
  facet_wrap(vars(Condition)) +
  labs(x = "RNA expression") +
  theme(legend.position = "none") +
  geom_vline(xintercept = 8000, col = "red", linetype = "dashed")

ggplot(rna, aes(x = Value, fill = Condition)) +
  geom_histogram() +
  xlim(c(0, 8000)) +
  ylim(c(0, 3000)) +
  facet_wrap(vars(Condition)) +
  labs(x = "RNA expression", title = "Zoomed distribution") +
  theme(legend.position = "none")








