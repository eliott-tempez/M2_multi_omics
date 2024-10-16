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

# See if there are nas
sum(is.na(data$RNA))
sum(is.na(data$protein))
sum(is.na(data$cyto))

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

# Cytokine expression distribution
ggplot(cyto, aes(y = Value)) +
  geom_boxplot() +
  labs(title = "Cytokine expression distribution (boxplot)",
       y = "Expression level")

ggplot(cyto, aes(x = Value)) +
  geom_histogram() +
  labs(title = "Cytokine expression distribution (histogram)",
    x = "Expression level")

# Expression of each molecule
ggplot(cyto, aes(x = Cytokine, y = Value)) +
  geom_boxplot() +
  labs(title = "Individual Cytokine expression", 
       x = "Cytokine", 
       y = "Expression level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Cytokine expression distribution, by condition
ggplot(cyto, aes(x = Condition, y = Value, col = Condition)) +
  geom_boxplot() +
  labs(title = "Cytokine expression distribution by condition", 
       x = "Time of prelevement", 
       y = "Expression level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggplot(cyto, aes(x = Value, fill = Condition)) +
  geom_histogram() +
  facet_wrap(vars(Condition)) +
  labs(title = "Cytokine expression distribution by condition",
    x = "Expression level") +
  theme(legend.position = "none")

# Expression of each molecule, by condition
ggplot(cyto, aes(x = Cytokine, y = Value, col = Condition)) +
  geom_boxplot() +
  labs(title = "Individual Cytokine Expression by time of sampling", 
       x = "Cytokine", 
       y = "Expression level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


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

# Protein expression distribution
ggplot(prot, aes(y = Value)) +
  geom_boxplot() +
  labs(title = "Protein expression distribution (boxplot)",
       y = "Expression level")

ggplot(prot, aes(x = Value)) +
  geom_histogram() +
  labs(title = "Protein expression distribution (histogram)",
       x = "Expression level")

# Expression of each molecule
ggplot(prot, aes(x = Protein, y = Value)) +
  geom_boxplot() +
  labs(title = "Individual Protein expression", 
       x = "Protein", 
       y = "Expression level") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Protein expression distribution, by condition
ggplot(prot, aes(x = Condition, y = Value, col = Condition)) +
  geom_boxplot() +
  labs(title = "Protein expression distribution by condition", 
       x = "Time of prelevement", 
       y = "Expression level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggplot(prot, aes(x = Value, fill = Condition)) +
  geom_histogram() +
  facet_wrap(vars(Condition)) +
  labs(title = "Protein expression distribution by condition",
       x = "Expression level") +
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

# RNA expression distribution
ggplot(rna, aes(y = Value)) +
  geom_boxplot() +
  labs(title = "RNA expression distribution (boxplot)",
       y = "Expression level")

ggplot(rna, aes(x = Value)) +
  geom_histogram() +
  labs(title = "RNA expression distribution (histogram)",
       x = "Expression level") +
  geom_vline(xintercept = 8000, col = "red", linetype = "dashed")

ggplot(rna, aes(x = Value)) +
  geom_histogram() +
  labs(title = "RNA expression distribution (zoomed histogram)",
       x = "Expression level") +
  xlim(c(0, 8000)) +
  ylim(c(0, 3000))

# RNA expression distribution, by condition
ggplot(rna, aes(x = Condition, y = Value, col = Condition)) +
  geom_boxplot() +
  labs(title = "RNA expression distribution by condition", 
       x = "Time of prelevement", 
       y = "Expression level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggplot(rna, aes(x = Value, fill = Condition)) +
  geom_histogram() +
  facet_wrap(vars(Condition)) +
  labs(title = "RNA expression distribution by condition",
       x = "Expression level") +
  theme(legend.position = "none") +
  geom_vline(xintercept = 8000, col = "red", linetype = "dashed")

ggplot(rna, aes(x = Value, fill = Condition)) +
  geom_histogram() +
  facet_wrap(vars(Condition)) +
  labs(title = "RNA expression distribution by condition (zoomed)",
       x = "Expression level") +
  theme(legend.position = "none") +
  xlim(c(0, 8000)) +
  ylim(c(0, 3000))






