## this script compares Chat+ vs Chat- nuclei basted on GSE161621 BLUM et al
## data with author annotations

## load libraries
library("dplyr")
library("readr")
library("stringr")
library("tidyr")
library("future")
library("Seurat")
library("anndata")
library("Matrix")

## seurat parameters
mito <- "^mt-" ## "^mt-" if mouse, "^MT-" if human
up_feat <- 6000
low_feat <- 200
percent_mito <- 20
low_count <- 400
var_feat <- 2000
NPC <- 30 ## number of PCs of first PCAs
k_para <- 50 ## number of k paramaters in FindNeighbors default is 30
sig_dims <- 12 ## number of PCs of final PCA

## parallelize workflow
future::plan("multicore", workers = 72) # uses 48 CPU
options(future.globals.maxSize = 5000 * 1024^2) ## 5GB per worker

###################### reconstruct original whole Blum data ####################
# Path to your H5AD file
fp <- "/staging/leuven/stg_00130/Projects/Diana_Piol_Project/author_annotated_MN/"
h5ad_file <- paste0(fp, "h5d_dataset/allexpsctl.h5ad")

# Read the H5AD file
adata <- anndata::read_h5ad(h5ad_file)

# Extract the matrix
mat <- adata$X

### read cholinergic only
readRDS(file = paste0(fp, "lepichon_gitler_integrated.RDS")
        ) -> lepichon_gitler_integrated

## filter to only contain Gitler data and relevant columns
lepichon_gitler_integrated@meta.data %>%
  tibble::rownames_to_column(var = "cell_id") %>%
  dplyr::filter(group == "Gitler") %>%
  dplyr::select(cell_id, cell_class) %>%
  as.data.frame() -> gitler_df

## read CSV Blum metadata table and left join with cell subtype annotations
readr::read_csv(file = paste0(
  fp, "Blum_all.exps.annotations.csv")) %>%
  dplyr::rename(cell_id = names(.)[1]) %>%
  dplyr::left_join(gitler_df, by = "cell_id") %>%
  as.data.frame() -> Blum_all_exps_annotations

## load processed sparse matrix into Seurat
Seurat::CreateSeuratObject(counts = t(mat),
                           min.cells = 1,
                           min.features = 1) -> obj

## process metadata
obj@meta.data %>%
  tibble::rownames_to_column(var = "cell_id") %>%
  dplyr::select(cell_id) %>%
  dplyr::left_join(Blum_all_exps_annotations, by = "cell_id") %>%
  dplyr::mutate(
    is_neuron = dplyr::case_when(
      cell_type == "Cholinergic neurons" ~ "neuron",
      cell_type == "Excitatory neurons" ~ "neuron",
      cell_type == "Inhibitory neurons" ~ "neuron",
      TRUE ~ "other"),
    chat_status = dplyr::case_when(
      cell_class == "Skeletal MNs" ~ "Chatpos",
      is_neuron == "neuron" ~ "Chatneg",
      TRUE ~ "Chatneg")) %>%
  tibble::column_to_rownames(var = "cell_id") %>%
  as.data.frame() -> obj@meta.data

## subset only neurons
obj2 <- subset(seu, subset = is_neuron == "neuron")

## set ident
Idents(object = obj2) <- "chat_status"
levels(obj2) ## now our DE is based on Chatpos vs Chatneg groups

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

######################### run FindAllMarkers to get DEGs #######################
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
saveRDS(object = smg, file = paste0(fp, "GSE161621_Blum_author_ann_DE_final.rds"))