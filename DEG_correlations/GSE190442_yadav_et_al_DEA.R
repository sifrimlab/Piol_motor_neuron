## this script takes the data from GSE190442 Yadav, et al. 2021, and compares
## Chat+ vs Chat- groups to produce a list of DE genes

## load libraries
library("dplyr")
library("ggplot2")
library("readr")
library("Seurat")
library("future")

## parallelize workflow
future::plan("multicore", workers = 72) # uses 48 CPU
options(future.globals.maxSize = 5000 * 1024^2) ## 5GB per worker

## features
mito <- "^mt-" ## "^mt-" if mouse, "^MT-" if human
up_feat <- 6000
low_feat <- 200
percent_mito <- 20
low_count <- 400
var_feat <- 2000
NPC <- 50 ## number of PCs of first PCA
k_para <- 50 ## number of k paramaters in FindNeighbors default is 30
sig_dims <- 50 ## number of PCs of final PCA

## file path to
fp <- "/staging/leuven/stg_00130/Projects/Diana_Piol_Project/RDS/"

## read metadata
readr::read_csv(file = paste0(fp, "GSE190442_aggregated_metadata_postqc.csv")) %>%
  dplyr::rename(cell_id = names(.)[1]) %>%
  dplyr::mutate(cell_id = gsub("-", ".", cell_id),
                motor_neuron = dplyr::case_when(
                  subtype_annotation == "Motoneurons" ~ "chat_pos",
                  subtype_annotation != "Motoneurons" ~ "chat_neg")) %>%
  dplyr::select(-c(nCount_RNA, nFeature_RNA)) %>%
  as.data.frame() -> meta_data

## read post-QC counts
readr::read_csv(file = paste0(fp, "GSE190442_aggregated_counts_postqc.csv")) %>%
  dplyr::rename(cell_id = names(.)[1]) %>%
  tibble::column_to_rownames(var = "cell_id") -> df1

## create Seurat object
obj <- Seurat::CreateSeuratObject(counts = as.matrix(df1),
                                  min.cells = 1,
                                  min.features = 1)

## coerce metadata
obj@meta.data %>%
  tibble::rownames_to_column(var = "cell_id") %>%
  dplyr::left_join(meta_data, by = "cell_id") %>%
  as.data.frame() -> seu_meta
rownames(seu_meta) <- seu_meta$cell_id
obj@meta.data <- seu_meta

## subset so data contains only neurons
obj2 <- subset(obj, top_level_annotation == "Neurons")

## do the usual Seurat workflow stuff
obj2 %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = var_feat) %>%
  ScaleData() %>%
  RunPCA(npcs = NPC, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:sig_dims, k.param = k_para,
                compute.SNN = TRUE, nn.method = "rann", nn.eps = 0) %>%
  FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) %>%
  RunUMAP(reduction = "pca", dims = 1:sig_dims) -> seu

## run FindAllMarkers
seu <- SetIdent(object = seu, value = "motor_neuron")
levels(seu)
seu %>%
  FindAllMarkers(test.use = "wilcox") %>%
  as.data.frame() %>%
  dplyr::arrange(p_val_adj) %>%
  dplyr::select(gene, cluster, everything()) %>%
  dplyr::rename(group = cluster) %>%
  dplyr::filter(group == "chat_pos") %>%
  as.data.frame() -> smg
rownames(smg) <- NULL

#@ save DE list
saveRDS(object = smg, file = paste0(fp, "GSE190442_DE.rds"))
