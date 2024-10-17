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




#############################################################################
#################################   DIABLO   ################################
#############################################################################

##### Genes and proteins #####
# set a list of all the X dataframes
data <-list(RNA = rna_matrix, 
            protein = prot_matrix)
Y <- condition_group

# Set design matrix
design <- matrix(0.1, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) <- 0 # set diagonal to 0s

# Form basic DIABLO model
basic.diablo.model = block.splsda(X = data, Y = Y, ncomp = 5, design = design) 

# run component number tuning with repeated CV
perf.diablo = perf(basic.diablo.model, validation = 'Mfold', 
                   folds = 10, nrepeat = 10) 

png(paste0(output_fold, 'diablo/elbow_plot_2.png'))
plot(perf.diablo, col = colors[1:3]) # plot output of tuning
dev.off()

# set the optimal ncomp value
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] 
# show the optimal choice for ncomp for each dist metric
perf.diablo$choice.ncomp$WeightedVote 

# set grid of values for each component to test
test.keepX <- list(
  RNA = c(5, 10, 20, 30, 40, 50),
  protein = c(5, 10, 20, 30, 40, 50))

# run the feature selection tuning
tune.TCGA <- tune.block.splsda(X = data, Y = Y, ncomp = ncomp, 
                              test.keepX = test.keepX, design = design,
                              validation = 'Mfold', folds = 10, nrepeat = 1,
                              dist = "centroids.dist")
list.keepX = tune.TCGA$choice.keepX # set the optimal values of features to retain

# set the optimised DIABLO model
final.diablo.model = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                                  keepX = list.keepX, design = design)

# Plot model
png(paste0(output_fold, 'diablo/plot_2.png'))
plotDiablo(final.diablo.model, ncomp = 1, col.per.group = colors)
dev.off()



##### Genes, proteins and cytokines #####
# set a list of all the X dataframes
data <-list(RNA = rna_matrix, 
            protein = prot_matrix,
            cyto = cyto_matrix)
Y <- condition_group

# Set design matrix
design <- matrix(0.1, ncol = length(data), nrow = length(data), 
                 dimnames = list(names(data), names(data)))
diag(design) <- 0 # set diagonal to 0s

# Form basic DIABLO model
basic.diablo.model = block.splsda(X = data, Y = Y, ncomp = 5, design = design) 

# run component number tuning with repeated CV
perf.diablo = perf(basic.diablo.model, validation = 'Mfold', 
                   folds = 10, nrepeat = 10) 

png(paste0(output_fold, 'diablo/elbow_plot_3.png'))
plot(perf.diablo, col = colors[1:3]) # plot output of tuning
dev.off()

# set the optimal ncomp value
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"] 
# show the optimal choice for ncomp for each dist metric
perf.diablo$choice.ncomp$WeightedVote 

# set grid of values for each component to test
test.keepX <- list(
  RNA = c(5, 10, 20, 30, 40, 50),
  protein = c(5, 10, 20, 30, 40, 50),
  cyto = c(5, 10, 15, 20))

# run the feature selection tuning
tune.TCGA <- tune.block.splsda(X = data, Y = Y, ncomp = ncomp, 
                               test.keepX = test.keepX, design = design,
                               validation = 'Mfold', folds = 10, nrepeat = 1,
                               dist = "centroids.dist")
list.keepX = tune.TCGA$choice.keepX # set the optimal values of features to retain

# set the optimised DIABLO model
final.diablo.model = block.splsda(X = data, Y = Y, ncomp = ncomp, 
                                  keepX = list.keepX, design = design)

# Plot model
png(paste0(output_fold, 'diablo/plot_3.png'),
           width = 3000, height = 3000, units = "px", pointsize = 72)
plotDiablo(final.diablo.model, ncomp = 1, col.per.group = colors)
dev.off()



##### Circo and network plot #####
png(paste0(output_fold, 'diablo/circo_plot.png'),
    width = 3000, height = 3000, units = "px", pointsize = 72)
circosPlot(final.diablo.model, cutoff = 0.7, line = TRUE,
           color.blocks= colors[1:3],
           color.cor = c("chocolate3","grey20"), size.labels = 1.5)
dev.off()

png(paste0(output_fold, 'diablo/network_plot.png'),
    width = 3000, height = 3000, units = "px", pointsize = 50)
network(final.diablo.model, 
        blocks = c(1,2,3),
        color.node = colors[2:4], 
        cutoff = 0.7, 
        cex.node = 0.4)
dev.off()







