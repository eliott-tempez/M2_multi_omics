library(tidyr)
library(ggplot2)
library(tidyverse)
library(mixOmics)

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
# See if there are 0 values
sum(data$cyto == 0)
sum(data$protein == 0)
sum(data$RNA == 0)



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


# Update data to get rid of values under the threshold
cyto_to_del <- which(cv_cyto$coef_var < 0.4)
cyto_to_del <- colnames(data$cyto)[cyto_to_del]
cyto <- cyto %>%
  filter(!Cytokine %in% cyto_to_del)
length(unique(cyto$Cytokine))
# On passe de 62 à 52 cytokines



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
  geom_vline(xintercept = 0.2, col = "red", linetype = "dashed")
ggsave(paste0(output_fold, 'variation_coef/prot_cv_lim.png'))


# Update data to get rid of values under the threshold
prot_to_del <- which(cv_prot$coef_var < 0.2)
prot_to_del <- colnames(data$prot)[prot_to_del]
prot <- prot %>%
  filter(!Protein %in% prot_to_del)
length(unique(prot$Protein))
# On passe de 2329 à 1610 protéines



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
  geom_histogram(bins = 200) +
  labs(title = "RNA expression distribution (histogram)",
       x = "Expression level") +
  geom_vline(xintercept = 1000, col = "red", linetype = "dashed") +
  geom_vline(xintercept = -1000, col = "red", linetype = "dashed")
ggsave(paste0(output_fold, 'initial_visualisation/rna_hist.png'))

ggplot(rna, aes(x = Value)) +
  geom_histogram() +
  labs(title = "RNA expression distribution (zoomed histogram)",
       x = "Expression level") +
  xlim(c(-1000, 1000))
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
  geom_histogram(bins = 200) +
  facet_wrap(vars(Condition)) +
  labs(title = "RNA expression distribution by condition",
       x = "Expression level") +
  theme(legend.position = "none") +
  geom_vline(xintercept = -1000, col = "red", linetype = "dashed") +
  geom_vline(xintercept = 1000, col = "red", linetype = "dashed")
ggsave(paste0(output_fold, 'initial_visualisation/rna_hist_condition.png'))

ggplot(rna, aes(x = Value, fill = Condition)) +
  geom_histogram() +
  facet_wrap(vars(Condition)) +
  labs(title = "RNA expression distribution by condition (zoomed)",
       x = "Expression level") +
  theme(legend.position = "none") +
  xlim(c(-1000, 1000))
ggsave(paste0(output_fold, 'initial_visualisation/rna_hist_zoom_condition.png'))


# Calculate variation coef
cv_rna <- rna %>%
  group_by(RNA) %>%
  summarise(
    mean_value = mean(Value),
    sd_value = sd(Value))
# Plot variation coef distribution
ggplot(cv_rna, aes(x = sd_value)) +
  geom_histogram(binwidth = 10) +
  labs(title = "RNA standard deviation distribution (zoomed)",
       x = "Standard deviation") +
  xlim(c(0, 3500)) +
  ylim(c(0, 400))
ggsave(paste0(output_fold, 'variation_coef/rna_cv.png'))

ggplot(cv_rna, aes(x = sd_value)) +
  geom_histogram(binwidth = 10) +
  labs(title = "RNA standard deviation distribution (zoomed)",
       x = "Standard deviation") +
  xlim(c(0, 3500)) +
  ylim(c(0, 400)) +
  geom_vline(xintercept = 80, col = "red", linetype = "dashed")
ggsave(paste0(output_fold, 'variation_coef/rna_cv_lim.png'))


# Update data to get rid of values under the threshold
rna_to_del <- which(cv_rna$sd_value < 80)
rna_to_del <- colnames(data$RNA)[rna_to_del]
rna <- rna %>%
  filter(!RNA %in% rna_to_del)
length(unique(rna$RNA))
# On passe de 23855 à 1853 arn



##### Scale data #####
cyto_scaled <- cyto %>%
  group_by(Cytokine) %>%
  mutate(Value = scale(Value, center = TRUE, scale = TRUE)) %>%
  ungroup()
prot_scaled <- prot %>%
  group_by(Protein) %>%
  mutate(Value = scale(Value, center = TRUE, scale = TRUE)) %>%
  ungroup()
rna_scaled <- rna %>%
  group_by(RNA) %>%
  mutate(Value = scale(Value, center = TRUE, scale = TRUE)) %>%
  ungroup()


# Expression distribution, by condition
ggplot(cyto_scaled, aes(x = Condition, y = Value, col = Condition)) +
  geom_boxplot() +
  labs(title = "Scaled Cytokin expression distribution", 
       x = "Time of prelevement", 
       y = "Expression level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave(paste0(output_fold, 'scaled_visualisation/cyto_boxplot_condition.png'))

ggplot(prot_scaled, aes(x = Condition, y = Value, col = Condition)) +
  geom_boxplot() +
  labs(title = "Scaled protein expression distribution", 
       x = "Time of prelevement", 
       y = "Expression level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave(paste0(output_fold, 'scaled_visualisation/prot_boxplot_condition.png'))

ggplot(rna_scaled, aes(x = Condition, y = Value, col = Condition)) +
  geom_boxplot() +
  labs(title = "Scaled RNA expression distribution", 
       x = "Time of prelevement", 
       y = "Expression level") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
ggsave(paste0(output_fold, 'scaled_visualisation/rna_boxplot_condition.png'))


##### Save processed data #####
processed_data <- list("rna" = rna_scaled, "prot" = prot_scaled, "cyto" = cyto_scaled)
save(processed_data, file = './data/processed_data.RData')



# Most varying molecules
cv_cyto %>%
  arrange(desc(coef_var)) %>%
  head(3)

cv_prot %>%
  arrange(desc(coef_var)) %>%
  head(3)

cv_rna %>%
  arrange(desc(sd_value)) %>%
  head(3)



#############################################################################
##################################   ACP   ##################################
#############################################################################
##### Cytokines #####
cyto_wide <- cyto_scaled %>%
  dplyr::select(Participant, Condition, Cytokine, Value) %>%
  tidyr::pivot_wider(names_from = Cytokine, values_from = Value)
# Get rid of non numeric variables
cyto_matrix <- as.matrix(dplyr::select(cyto_wide, -Participant, -Condition))
condition_group <- cyto_wide$Condition

# Elbow plot
pca_cyto <- pca(cyto_matrix, ncomp = 10)
plot(pca_cyto)
ggsave(paste0(output_fold, 'pca/elbow_plot_cyto.png'))
# On choisit 3 composantes
pca_cyto <- pca(cyto_matrix, ncomp = 3)

# Plot PCA
plotIndiv(pca_cyto, comp=c(1, 2), group = condition_group)
ggsave(paste0(output_fold, 'pca/pca_cyto_indiv.png'))
plotVar(pca_cyto)
ggsave(paste0(output_fold, 'pca/pca_cyto_var.png'))


##### Proteins #####
prot_wide <- prot_scaled %>%
  dplyr::select(Participant, Condition, Protein, Value) %>%
  tidyr::pivot_wider(names_from = Protein, values_from = Value)
# Get rid of non numeric variables
prot_matrix <- as.matrix(dplyr::select(prot_wide, -Participant, -Condition))

# Elbow plot
pca_prot <- pca(prot_matrix, ncomp = 10)
plot(pca_prot)
ggsave(paste0(output_fold, 'pca/elbow_plot_prot.png'))
# On choisit 3 composantes
pca_prot <- pca(prot_matrix, ncomp = 3)

# PCA plot
plotIndiv(pca_prot, group = condition_group)
ggsave(paste0(output_fold, 'pca/pca_prot_indiv.png'))
plotVar(pca_prot)
ggsave(paste0(output_fold, 'pca/pca_prot_indiv.png'))


##### RNA #####
rna_wide <- rna_scaled %>%
  dplyr::select(Participant, Condition, RNA, Value) %>%
  tidyr::pivot_wider(names_from = RNA, values_from = Value)
# Get rid of non numeric variables
rna_matrix <- as.matrix(dplyr::select(rna_wide, -Participant, -Condition))

# Elbow plot
pca_rna <- pca(rna_matrix, ncomp = 10)
plot(pca_rna)
ggsave(paste0(output_fold, 'pca/elbow_plot_rna.png'))
# On choisit 2 composantes
pca_rna <- pca(rna_matrix, ncomp = 2)

# PCA plot
plotIndiv(pca_rna, group = condition_group)
ggsave(paste0(output_fold, 'pca/pca_prot_indiv.png'))
plotVar(pca_rna)
ggsave(paste0(output_fold, 'pca/pca_prot_indiv.png'))

# Contributing variables
# Get the loadings for the first two principal components
loadings_pca <- pca_rna$loadings$X
head(loadings_pca[, 1])
head(loadings_pca[, 2])












