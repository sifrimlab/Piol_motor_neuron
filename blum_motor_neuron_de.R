## new analysis of GSE161621 BLUM et al data with author annotations

## load libraries
library("dplyr")
library("readr")
library("stringr")
library("tidyr")
library("future")
library("Seurat")
library("anndata")
library("Matrix")

## ensure Seurat assay is V3 and not V5
# options(Seurat.object.assay.version = "v3")
###################### reconstruct original whole Blum data ####################
# Path to H5AD file
# fp <- "/staging/leuven/stg_00130/Projects/Diana_Piol_Project/author_annotated_MN/"
fp <- "~/Documents/tmp/dacruz/202403_piol/data/author_annotated_MN/"
h5ad_file <- paste0(fp, "h5d_dataset/allexpsctl.h5ad")

# Read the H5AD file
adata <- anndata::read_h5ad(h5ad_file)

# Extract the matrix
mat <- adata$X

### read cholinergic only
readRDS(file = paste0(fp, "lepichon_gitler_integrated.RDS")
        ) -> lepichon_gitler_integrated
# View(lepichon_gitler_integrated@meta.data)

## filter to only contain Gitler data and relevant columns
lepichon_gitler_integrated@meta.data %>%
  tibble::rownames_to_column(var = "cell_id") %>%
  dplyr::filter(group == "Gitler") %>%
  dplyr::select(cell_id, cell_class) %>%
  as.data.frame() -> blum_df

## read CSV Blum metadata table and left join with cell subtype annotations
readr::read_csv(file = paste0(
  fp, "Blum_all.exps.annotations.csv")) %>%
  dplyr::rename(cell_id = names(.)[1]) %>%
  dplyr::left_join(blum_df, by = "cell_id") %>%
  as.data.frame() -> Blum_all_exps_annotations
# View(Blum_all_exps_annotations)

## load processed sparse matrix into Seurat
Seurat::CreateSeuratObject(counts = t(mat),
                           assay = "RNA",
                           min.cells = 1,
                           min.features = 1) -> seu
# View(seu@meta.data)

## process metadata
seu@meta.data %>%
  tibble::rownames_to_column(var = "cell_id") %>%
  dplyr::select(cell_id) %>%
  dplyr::left_join(Blum_all_exps_annotations, by = "cell_id") %>%
  dplyr::mutate(
    is_neuron = dplyr::case_when(
      cell_type == "Cholinergic neurons" ~ "neuron",
      cell_type == "Excitatory neurons" ~ "neuron",
      cell_type == "Inhibitory neurons" ~ "neuron",
      TRUE ~ "other"), ## differentiate between neurons and non-neurons
    chat_status = dplyr::case_when(
      cell_class == "Skeletal MNs" ~ "Chatpos",
      is_neuron == "neuron" ~ "Chatneg",
      TRUE ~ "Chatneg")) %>%
  tibble::column_to_rownames(var = "cell_id") %>%
  as.data.frame() -> seu@meta.data

# ## subset only neurons
# neuron_seu <- subset(seu, subset = is_neuron == "neuron")
# 
# ## set ident
# Idents(object = neuron_seu) <- "chat_status"
# levels(neuron_seu) ## now our DE is based on Chatpos vs Chatneg groups
# # table(neuron_seu@meta.data$chat_status)
# # table(neuron_seu@meta.data$cell_type)
# # table(neuron_seu@meta.data$cell_class)
# saveRDS(object = neuron_seu, file = paste0(fp, "blum_gitler_neurons_correct-2.rds"))

###############################################################################

## parallelize workflow
future::plan("multicore", workers = 72) # uses 48 CPU
options(future.globals.maxSize = 5000 * 1024^2) ## 5GB per worker

## correct
# Chatneg Chatpos
# 10284    3305

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
saveRDS(object = smg, file = paste0(fp, "GSE161621_Blum_author_ann_DE_correct.rds"))

############################## inclusion analysis ##############################
library("VennDiagram")
### Blum new + old DE inclusion analysis
readRDS("~/Documents/tmp/dacruz/202403_piol/data/new_DE/GSE161621_Blum_author_ann_DE_final.rds") %>%
  dplyr::rename(genes = names(.)[1]) %>%
  as.data.frame() -> GSE161621_Blum_author_ann_DE_final
# View(GSE161621_Blum_author_ann_DE_final)

## write new DE results to Excel
local_path <- "~/Documents/tmp/dacruz/202403_piol/data/new_DE/"
writexl::write_xlsx(x = GSE161621_Blum_author_ann_DE_final,
                    path = paste0(
                      local_path, "GSE161621_Blum_author_ann_Chatpos_v_Chat_neg_DE_final_res.xlsx"))