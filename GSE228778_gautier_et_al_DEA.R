## this script takes the data from GSE228778 Gautier, et al. 2023, and compares
## Chat+ vs Chat- groups to produce a list of DE genes

## load libraries
library("dplyr")
library("tibble")
library("ggplot2")
library("Seurat")
library("future")
library("writexl")

## parallelize workflow
future::plan("multicore", workers = 72) # uses 48 CPU
options(future.globals.maxSize = 5000 * 1024^2) ## 5GB per worker

## file path
fp <- "/staging/leuven/stg_00130/Projects/Diana_Piol_Project/GSE228778/"

## NOTE: this dataset is already a Seurat object with normalized and scaled
## counts, so we can load the complete RDS and proceed directly to perform DEG
## analysis
seu <- readRDS(file = paste0(fp, "GSE228778_gautier_neurons.rds"))

## set levels
seu <- SetIdent(object = seu, value = "motor_neuron")
# levels(seu)

## run FindAllMarkers
seu %>%
  FindAllMarkers(test.use = "wilcox") %>%
  as.data.frame() %>%
  dplyr::arrange(p_val_adj) %>%
  tibble::rownames_to_column(var = "rowname") %>%
  dplyr::select(gene, cluster, everything()) %>%
  dplyr::select(-rowname) %>%
  dplyr::rename(group = cluster) %>%
  dplyr::filter(group == "Motor Neurons") %>%
  as.data.frame() -> deg_list

## save motor neuron genes
saveRDS(object = deg_list, file = paste0(fp, "GSE228778_motoneuron_DEGs.rds"))