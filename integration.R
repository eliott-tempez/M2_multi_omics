library(tidyr)
library(ggplot2)
library(tidyverse)
library(mixOmics)

# Set default image width
output_fold <- './output/plots/'

# Load data
load('./data/processed_data.RData')
cyto_scaled <- processed_data$cyto
prot_scaled <- processed_data$prot
rna_scaled <- processed_data$rna

#############################################################################
##################################   PLS   ##################################
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


##### PLS #####
colors <- c("1" = "#3d1619", "2" = "#69150f", "3" = "#a44200", "4" = "#d58937")

# With 10 composants
pls <- pls(X = rna_matrix, Y = prot_matrix, ncomp = 10)
# Calcul de la proportion de covariance expliquée pour chaque composante
explained_var_X <- pls$prop_expl_var$X
explained_var_Y <- pls$prop_expl_var$Y
# Moyenne de la proportion de covariance entre X et Y pour chaque composante
explained_var <- rowMeans(cbind(explained_var_X, explained_var_Y))

# Tracer la proportion de covariance capturée
png(paste0(output_fold, 'pls/elbow_plot.png'))
plot(explained_var, type = "b", pch = 19, col = "#3d1619", 
     xlab = "Nombre de composantes", ylab = "Proportion de Covariance Capturée",
     main = "Proportion de Covariance Capturée - PLS")
dev.off()
ggsave(paste0(output_fold, 'pls/elbow_plot.png'))
# On choisit 5 composantes
pls <- pls(X = rna_matrix, Y = prot_matrix, ncomp = 5)

# Plot PLS
plotIndiv(pls, group = condition_group, legend = TRUE, col.per.group = colors,
          title = "Graphe des Échantillons\nARN et Protéines")
ggsave(paste0(output_fold, 'pls/pls_indiv.png'))

plotVar(pls, col = c("#a44200", "#d58937"),
        blocks = c("RNA", "Protein"), legend = c("RNA", "Protein"))
ggsave(paste0(output_fold, 'pls/pls_var.png'))

plotArrow(pls, 
          title = "Arrow Plot - ARN et Protéines",
          legend = TRUE, legend.title = "Condition",
          col.per.group = colors, group = condition_group)
ggsave(paste0(output_fold, 'pls/arrow_plot.png'))





