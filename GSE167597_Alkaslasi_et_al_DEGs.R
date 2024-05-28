## this script compares Chat+ vs Chat- nuclei basted on GSE167597 Alkaslasi et al
## data with author annotations
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167597

## NOTE:: scaling and normalization of the counts is not necessary because we are
## working directly from the scaled and normalized author counts

## load libraries
library("dplyr")
library("readr")
library("stringr")
library("tidyr")
library("future")
library("Seurat")
library("anndata")
library("Matrix")

## parallelize workflow
future::plan("multicore", workers = 72) # uses 48 CPU
options(future.globals.maxSize = 5000 * 1024^2) ## 5GB per worker

## file path
fp <- "/staging/leuven/stg_00130/Projects/Diana_Piol_Project/author_annotated_MN/"

## read the RDS object for the full data
seu <- readRDS(file = paste0(fp, "AlkaslasiPiccusetal_allnuclei.RDS")) 

## change metadata to create Chatpos/Chatneg variables - everything is a neuron
## chat pos = Skeletal Motor Neurons
## chat neg = neurons

seu@meta.data %>% 
  dplyr::mutate(
    chat_status = dplyr::case_when(
      cholinergictypes == "Skeletal Motor Neurons" ~ "Chatpos",
      TRUE ~ "Chatneg")) %>% ## otherwise
  as.data.frame() -> seu@meta.data

## set ident
Idents(object = seu) <- "chat_status"
levels(seu) ## now our DE is based on Chatpos vs Chatneg groups
# table(seu@meta.data$chat_status)

## correct
# Chatneg Chatpos 
# 28856    5375 

## run FindAllMarkers to get DEGs
seu %>%
  FindAllMarkers(test.use = "wilcox") %>%
  as.data.frame() %>%
  dplyr::arrange(p_val_adj) %>%
  dplyr::select(gene, cluster, everything()) %>%
  dplyr::rename(group = cluster) %>%
  dplyr::filter(group == "Chatpos") %>%
  tibble::rownames_to_column(var = "row_names") %>%
  dplyr::select(-row_names) %>%
  as.data.frame() -> smg

#@ save DE list
saveRDS(object = smg, file = paste0(fp, "GSE167597_Alkaslasi_author_ann_v4_DE.rds"))