library(tidyr)
library(ggplot2)
library(tidyverse)

# Set default image width
output_fold <- './output/plots/'

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
# See if there are negative expression levels
sum(data$cyto < 0)
sum(data$protein < 0)
sum(data$RNA < 0)
paste("Number of negative cells:", length(which(data$RNA < 0, arr.ind = T)))

##### CYTOKINES #####
# Transform data
cyto <- data$cyto %>%
  # Add a new column splitting the rownames into "Participant" and "Condition"
  mutate(Participant = sub("_.*", "", rownames(data$cyto)),
         Condition = sub(".*_", "", rownames(data$cyto))) %>%
  # Reshape from wide to long format
  pivot_longer(cols = -c(Participant, Condition), # Pivot all columns except Participant and Condition
               names_to = "Cytokine",             # The names of the columns will go into "Cytokine"
               values_to = "Value")

# Cytokine expression distribution
ggplot(cyto, aes(y = Value)) +
  geom_boxplot() +
  labs(title = "Cytokine expression distribution (boxplot)",
       y = "Expression level")
ggsave(paste0(output_fold, 'initial_visualisation/cyto_boxplot.png'))

ggplot(cyto, aes(x = Value)) +
  geom_histogram() +
  labs(title = "Cytokine expression distribution (histogram)",
    x = "Expression level")
ggsave(paste0(output_fold, 'initial_visualisation/cyto_hist.png'))

# Expression of each molecule
ggplot(cyto, aes(x = Cytokine, y = Value)) +
  geom_boxplot() +
  labs(title = "Individual Cytokine expression", 
       x = "Cytokine", 
       y = "Expression level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(output_fold, 'initial_visualisation/cyto_boxplot_indiv.png'))

# Cytokine expression distribution, by condition
ggplot(cyto, aes(x = Condition, y = Value, col = Condition)) +
  geom_boxplot() +
  labs(title = "Cytokine expression distribution by condition", 
       x = "Time of prelevement", 
       y = "Expression level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave(paste0(output_fold, 'initial_visualisation/cyto_boxplot_condition.png'))

ggplot(cyto, aes(x = Value, fill = Condition)) +
  geom_histogram() +
  facet_wrap(vars(Condition)) +
  labs(title = "Cytokine expression distribution by condition",
    x = "Expression level") +
  theme(legend.position = "none")
ggsave(paste0(output_fold, 'initial_visualisation/cyto_hist_condition.png'))

# Expression of each molecule, by condition
ggplot(cyto, aes(x = Cytokine, y = Value, col = Condition)) +
  geom_boxplot() +
  labs(title = "Individual Cytokine Expression by time of sampling", 
       x = "Cytokine", 
       y = "Expression level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(output_fold, 'initial_visualisation/cyto_boxplot_indiv_condition.png'))


# Calculate variation coef
cv_cyto <- cyto %>%
  group_by(Cytokine) %>%
  summarise(
    mean_value = mean(Value),
    sd_value = sd(Value),
    coef_var = sd_value / mean_value
  )
# Plot variation coef distribution
ggplot(cv_cyto, aes(x = coef_var)) +
  geom_histogram() +
  labs(title = "Cytokine coef of variation distribution",
    x = "Expression coefficient of Variation (CV)")
ggsave(paste0(output_fold, 'variation_coef/cyto_cv.png'))

ggplot(cv_cyto, aes(x = coef_var)) +
  geom_histogram() +
  labs(title = "Cytokine coef of variation distribution",
       x = "Expression coefficient of Variation (CV)") +
  geom_vline(xintercept = 0.4, col = "red", linetype = "dashed")
ggsave(paste0(output_fold, 'variation_coef/cyto_cv_lim.png'))


##### PROTEINS #####
# Transform data
prot <- data$protein %>%
  # Add a new column splitting the rownames into "Participant" and "Condition"
  mutate(Participant = sub("_.*", "", rownames(data$cyto)),
         Condition = sub(".*_", "", rownames(data$cyto))) %>%
  # Reshape from wide to long format
  pivot_longer(cols = -c(Participant, Condition), # Pivot all columns except Participant and Condition
               names_to = "Protein",
               values_to = "Value")

# Protein expression distribution
ggplot(prot, aes(y = Value)) +
  geom_boxplot() +
  labs(title = "Protein expression distribution (boxplot)",
       y = "Expression level")
ggsave(paste0(output_fold, 'initial_visualisation/prot_boxplot.png'))

ggplot(prot, aes(x = Value)) +
  geom_histogram() +
  labs(title = "Protein expression distribution (histogram)",
       x = "Expression level")
ggsave(paste0(output_fold, 'initial_visualisation/prot_hist.png'))

# Expression of each molecule
ggplot(prot, aes(x = Protein, y = Value)) +
  geom_boxplot() +
  labs(title = "Individual Protein expression", 
       x = "Protein", 
       y = "Expression level") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
ggsave(paste0(output_fold, 'initial_visualisation/prot_boxplot_indiv.png'))

# Protein expression distribution, by condition
ggplot(prot, aes(x = Condition, y = Value, col = Condition)) +
  geom_boxplot() +
  labs(title = "Protein expression distribution by condition", 
       x = "Time of prelevement", 
       y = "Expression level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave(paste0(output_fold, 'initial_visualisation/prot_boxplot_condition.png'))

ggplot(prot, aes(x = Value, fill = Condition)) +
  geom_histogram() +
  facet_wrap(vars(Condition)) +
  labs(title = "Protein expression distribution by condition",
       x = "Expression level") +
  theme(legend.position = "none")
ggsave(paste0(output_fold, 'initial_visualisation/prot_hist_condition.png'))


# Calculate variation coef
cv_prot <- prot %>%
  group_by(Protein) %>%
  summarise(
    mean_value = mean(Value),
    sd_value = sd(Value),
    coef_var = sd_value / mean_value
  )
# Plot variation coef distribution
ggplot(cv_prot, aes(x = coef_var)) +
  geom_histogram() +
  labs(title = "Protein coef of variation distribution",
       x = "Expression coefficient of Variation (CV)")
ggsave(paste0(output_fold, 'variation_coef/prot_cv.png'))

# Plot variation coef distribution
ggplot(cv_prot, aes(x = coef_var)) +
  geom_histogram() +
  labs(title = "Protein coef of variation distribution",
       x = "Expression coefficient of Variation (CV)") +
  geom_vline(xintercept = 0.125, col = "red", linetype = "dashed")
ggsave(paste0(output_fold, 'variation_coef/prot_cv_lim.png'))


##### RNA #####
# Transform data
rna <- data$RNA %>%
  # Add a new column splitting the rownames into "Participant" and "Condition"
  mutate(Participant = sub("_.*", "", rownames(data$cyto)),
         Condition = sub(".*_", "", rownames(data$cyto))) %>%
  # Reshape from wide to long format
  pivot_longer(cols = -c(Participant, Condition), # Pivot all columns except Participant and Condition
               names_to = "RNA",
               values_to = "Value")

# RNA expression distribution
ggplot(rna, aes(y = Value)) +
  geom_boxplot() +
  labs(title = "RNA expression distribution (boxplot)",
       y = "Expression level")
ggsave(paste0(output_fold, 'initial_visualisation/rna_boxplot.png'))

ggplot(rna, aes(x = Value)) +
  geom_histogram() +
  labs(title = "RNA expression distribution (histogram)",
       x = "Expression level") +
  geom_vline(xintercept = 8000, col = "red", linetype = "dashed")
ggsave(paste0(output_fold, 'initial_visualisation/rna_hist.png'))

ggplot(rna, aes(x = Value)) +
  geom_histogram() +
  labs(title = "RNA expression distribution (zoomed histogram)",
       x = "Expression level") +
  xlim(c(0, 8000)) +
  ylim(c(0, 3000))
ggsave(paste0(output_fold, 'initial_visualisation/rna_hist_zoom.png'))

# RNA expression distribution, by condition
ggplot(rna, aes(x = Condition, y = Value, col = Condition)) +
  geom_boxplot() +
  labs(title = "RNA expression distribution by condition", 
       x = "Time of prelevement", 
       y = "Expression level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave(paste0(output_fold, 'initial_visualisation/rna_boxplot_condition.png'))

ggplot(rna, aes(x = Value, fill = Condition)) +
  geom_histogram() +
  facet_wrap(vars(Condition)) +
  labs(title = "RNA expression distribution by condition",
       x = "Expression level") +
  theme(legend.position = "none") +
  geom_vline(xintercept = 8000, col = "red", linetype = "dashed")
ggsave(paste0(output_fold, 'initial_visualisation/rna_hist_condition.png'))

ggplot(rna, aes(x = Value, fill = Condition)) +
  geom_histogram() +
  facet_wrap(vars(Condition)) +
  labs(title = "RNA expression distribution by condition (zoomed)",
       x = "Expression level") +
  theme(legend.position = "none") +
  xlim(c(0, 8000)) +
  ylim(c(0, 3000))
ggsave(paste0(output_fold, 'initial_visualisation/rna_hist_zoom_condition.png'))


# Calculate variation coef
cv_rna <- rna %>%
  group_by(RNA) %>%
  summarise(
    mean_value = mean(Value),
    sd_value = sd(Value),
    coef_var = abs(sd_value / mean_value)
  )
# Plot variation coef distribution
ggplot(cv_rna, aes(x = coef_var)) +
  geom_histogram() +
  labs(title = "RNA coef of variation distribution",
       x = "Expression coefficient of Variation (CV)")
ggsave(paste0(output_fold, 'variation_coef/rna_cv.png'))

ggplot(cv_rna, aes(x = coef_var)) +
  geom_histogram() +
  labs(title = "RNA coef of variation distribution",
       x = "Expression coefficient of Variation (CV)") +
  geom_vline(xintercept = 0.6, col = "red", linetype = "dashed")
ggsave(paste0(output_fold, 'variation_coef/rna_cv_lim.png'))





