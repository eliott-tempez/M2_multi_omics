library(tidyr)
library(ggplot2)
library(tidyverse)
library(mixOmics)

# Set default image width
output_fold <- './output/plots/'

# Load data
load('./data/processed_data.RData')
cyto <- processed_data$cyto
prot <- processed_data$prot
rna <- processed_data$rna