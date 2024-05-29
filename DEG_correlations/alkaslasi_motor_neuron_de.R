## new analysis of GSE167597 Alkaslasi data with new author annotations

## load libraries
library("dplyr")
library("readr")
library("stringr")
library("tidyr")
library("future")
library("Seurat")
library("readxl")

seu <- readRDS("~/Documents/tmp/dacruz/202403_piol/data/author_annotated_MN/AlkaslasiPiccusetal_allnuclei.RDS")
View(seu@meta.data)
table(seu@meta.data$cholinergictypes)

## read all markers
readxl::read_excel("~/Downloads/41467_2021_22691_MOESM4_ESM.xls") %>%
  dplyr::rename(row_names = names(.)[1]) %>%
  dplyr::select(-row_names) %>%
  as.data.frame() -> mol_DEGs
# View(mol_DEGs)
# unique(mol_DEGs$cluster) ## 0-37

## markers
readr::read_csv("~/Downloads/41467_2021_22691_MOESM5_ESM.csv") %>%
  dplyr::rename(row_names = names(.)[1]) %>%
  dplyr::select(-row_names) %>%
  dplyr::select(gene, everything()) %>%
  dplyr::arrange(p_val_adj) %>%
  as.data.frame() -> mol_DEGs_MN
View(mol_DEGs_MN)
unique(mol_DEGs_MN$cluster) ## 0-7

## parallelize workflow
future::plan("multicore", workers = 72) # uses 48 CPU
options(future.globals.maxSize = 5000 * 1024^2) ## 5GB per worker

# GSE167597
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167597
fp <- "/staging/leuven/stg_00130/Projects/Diana_Piol_Project/author_annotated_MN/"

## read the RDS object for the full data
readRDS(file = paste0(fp, "AlkaslasiPiccusetal_allnuclei.RDS")) -> seu

## various stats
# View(seu@meta.data)
# unique(seu@meta.data$types)
# table(seu@meta.data$types)

## everything is a neuron
## chat pos = Skeletal Motor Neurons
## chat neg = neurons

## change metadata to create Chatpos/Chatneg variables
seu@meta.data %>% 
  dplyr::mutate(
    chat_status = dplyr::case_when(
      cholinergictypes == "Skeletal Motor Neurons" ~ "Chatpos",
      TRUE ~ "Chatneg" ## all other neurons are Chatneg
      )) %>% 
  as.data.frame() -> seu@meta.data

## set ident
Idents(object = seu) <- "chat_status"
levels(seu) ## now our DE is based on Chatpos vs Chatneg groups
table(seu@meta.data$chat_status)

# gene_list <- c("Alcam", "Bnc2", "C1qtnf4", "Cpne4", "Dach2", "Erbb4", "Foxp1",
#                "Glis3", "Gpc3", "Gpr149", "Grm5", "Isl1", "Isl2", "Lhx3", "Mnx1",
#                "Mpped2", "Nrp2", "Piezo2", "Pitx2", "Plekhg1", "Rbfox3", "Reln",
#                "Rreb1", "Sema5a", "Stk32a", "Sst", "Sv2b", "Tox")
# head(AverageExpression(object = seu, features = gene_list)$RNA)

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