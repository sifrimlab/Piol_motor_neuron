## This script prepares count matrices of DAPI stained *spinal cord* RESOLVE
## ALS data and analyzes this data using Seurat for  samples

# Can you also do the same but with these comparisons:
# - control FUS vs FUS
# - control TDP vs TDP
# - control FUS & control TDP vs TDP
# - control FUS & control TDP vs FUS
# - control SOD vs SOD

## load libraries
library("dplyr")
library("readr")
library("writexl")
library("Seurat")
library("tibble")
library("ggplot2")
library("future")

## future
# future::plan("multicore", workers = 12) # uses 48 CPU
# options(future.globals.maxSize = 5000 * 1024^2) ## 5GB per worker

##################################### SC ######################################
################################# "SC DAPI_trans" #############################

## process Excel count matrices

## file path results
fp0 <- "~/Documents/tmp/dacruz/202307_piol/"

## file path input files
fp <- paste0("~/Documents/tmp/dacruz/202305_piol/count_matrices_correct/with_area/",
             "SC_Resolve_Images/")
lf <- list.files(path = fp, pattern = "*.csv")
dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    dplyr::filter(area > 2000, ## area QC filters
                  area < 15000) %>%
    dplyr::select(- c(area, tidyr::contains("FP"))) %>% ## remove area + FPs!
    tibble::column_to_rownames(var = "cell_label")})
dfs3 <- lapply(dfs2, function(i){t(i)})

## create list of seurat objects from the transformed count matrices
fls1 <- lapply(dfs3, function(i) CreateSeuratObject(counts = i,
                                                    min.cells = 1,
                                                    min.features = 1))
fls2 <- list()
for (i in 1:length(fls1)) {
  ## get experimental control condition and add to metadata
  fls1[[i]]@meta.data %>%
    dplyr::mutate(orig.ident = gsub("_.*", "", lf)[i]) -> fls2[[i]]
  fls1[[i]]@meta.data <- fls2[[i]]         ## add sample-wise metadata to object
  fls1[[i]]@meta.data$sample <- paste(i)   ## add sample number
  ## add x and y coords
  fls1[[i]]@meta.data$y <- as.numeric(
    gsub("row", "", gsub("_.*", "", rownames(fls2[[i]]))))
  fls1[[i]]@meta.data$x <- as.numeric(
    gsub("col", "", gsub(".*_", "", gsub("_[^_]+$", "", rownames(fls2[[i]])))))
}

## merge list of processed Excel count matrices into single Seurat object
merge(x = fls1[[1]], y = fls1[2:length(fls1)]) %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst") %>%
  ScaleData() %>%
  RunPCA(npcs = 10, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:10, compute.SNN = TRUE,
                nn.method = "rann", nn.eps = 0) %>%
  FindClusters(resolution = c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)) %>%
  RunUMAP(reduction = "pca", dims = 1:10) -> seu

## save prepared Seurat object
saveRDS(object = seu, file = paste0(fp0, "SC_filtered_seurat.rds"))
# 
# ## plots plots plots
# Idents(seu) <- "orig.ident"
# UMAPPlot(seu, reduction = "umap", pt.size = 0.9, label = FALSE) +
#   labs(title = paste0(
#     "UMAP of SC DAPI segmented objects by exp. condition")) +
#   theme(plot.title = element_text(hjust = 0.5)) -> p1
# 
# Idents(seu) <- "RNA_snn_res.0.1"
# UMAPPlot(seu, reduction = "umap", pt.size = 0.9, label = FALSE) +
#   labs(title = paste0(
#     "UMAP of SC DAPI at resolution 0.1")) +
#   theme(plot.title = element_text(hjust = 0.5)) -> p2
# 
# Idents(seu) <- "RNA_snn_res.0.2"
# UMAPPlot(seu, reduction = "umap", pt.size = 0.9, label = FALSE) +
#   labs(title = paste0(
#     "UMAP of SC DAPI at resolution 0.2")) +
#   theme(plot.title = element_text(hjust = 0.5)) -> p3
# 
# Idents(seu) <- "RNA_snn_res.0.3"
# UMAPPlot(seu, reduction = "umap", pt.size = 0.9, label = FALSE) +
#   labs(title = paste0(
#     "UMAP of SC DAPI at resolution 0.3")) +
#   theme(plot.title = element_text(hjust = 0.5)) -> p4
# 
# Idents(seu) <- "RNA_snn_res.0.4"
# UMAPPlot(seu, reduction = "umap", pt.size = 0.9, label = FALSE) +
#   labs(title = paste0(
#     "UMAP of SC DAPI at resolution 0.4")) +
#   theme(plot.title = element_text(hjust = 0.5)) -> p5
# 
# Idents(seu) <- "RNA_snn_res.0.5"
# UMAPPlot(seu, reduction = "umap", pt.size = 0.9, label = FALSE) +
#   labs(title = paste0(
#     "UMAP of SC DAPI at resolution 0.5")) +
#   theme(plot.title = element_text(hjust = 0.5)) -> p6
# 
# ## motor neuron marker genes
# list("Chat", # ChAT	general MN marker
#      "Cdh6",  # Cdh6	fast-firing MN marker
#      "Chodl", # Chodl	fast-firing MN marker
#      "Kcnd2", # Kcnt2	fast-firing MN marker
#      "Kcnq5", # Kcnq5	fast-firing MN marker
#      "Prkcb", # Prkcb	fast-firing MN marker
#      "Kcnt2", # Kcnd2	slow-firing MN marker
#      "Prkcd", # Prkcd	slow-firing MN marker
#      "Sv2a" # Sv2a	slow-firing MN marker
#      ) -> m1
# 
# DefaultAssay(seu) <- "RNA"
# p7 <- FeaturePlot(seu, features = m1, reduction = "umap", raster = FALSE)
# 
# Idents(seu) <- "RNA_snn_res.0.2"
# pl1 <- list(); pl2 <- list()
# for (i in 1:length(m1)){
# FeaturePlot(object = seu, features = m1[[i]], reduction = "umap",
#             cols = c("azure3", "darkblue"), pt.size = 0.7, raster = FALSE
#             ) -> pl1[[i]]
# VlnPlot(object = seu, features = m1[[i]], pt.size = 0) +
#   theme(axis.text.x = element_text(angle = 45, size = 10)) -> pl2[[i]]
# }
# 
# pdf(file = paste0(fp0, "new_spinal_cord_UMAPs.pdf"))
# p1; p2; p3; p4; p5; p6; p7;
# for (i in 1:length(m1)){
# print(pl1[[i]])
# print(pl2[[i]])
# }
# dev.off()
# 
# ## QC plots
# seu@meta.data %>%
#   as.data.frame() %>%
#   ggplot(aes(x = nCount_RNA)) +
#   geom_density() +
#   ggtitle("Distribution of total counts per segmented cell in SC DAPI") -> q1
# seu@meta.data %>%
#   as.data.frame() %>%
#   ggplot(aes(x = nFeature_RNA)) +
#   geom_density() +
#   ggtitle("Distribution of total genes per segmented cell in SC DAPI") -> q2
# seu@meta.data %>%
#   as.data.frame() %>%
#   ggplot(aes(x = x, y = y, color = nCount_RNA)) +
#   geom_point(size = 0.1) +
#   ggtitle("Distribution of nCount per segmented object in SC DAPI") +
#   scale_colour_gradient(low = "azure3", high = "red" ) -> q3
# seu@meta.data %>%
#   as.data.frame() %>%
#   ggplot(aes(x = x, y = y, color = nFeature_RNA)) +
#   geom_point(size = 0.1) +
#   ggtitle("Distribution of nFeature per segmented object in SC DAPI") +
#   scale_colour_gradient(low = "azure3", high = "red" ) -> q4
# 
# pdf(file = paste0(fp0, "new_spinal_cord_qc_plots.pdf"))
# q1; q2; q3; q4
# dev.off()

################################ spatial plots #################################
# seu@meta.data %>%
#   dplyr::mutate(sample_ident = paste0(orig.ident, "_", sample)) %>%
#   as.data.frame() -> seu@meta.data
# 
# gpl <- list()
# for (i in 1:length(unique(seu@meta.data$sample_ident))) {
# seu@meta.data %>%
#   dplyr::filter(sample_ident == unique(seu@meta.data$sample_ident)[i]) %>%
#   ggplot(aes(x = x, y = y, color = RNA_snn_res.0.2)) +
#   geom_point(size = 1.0, alpha = 1) +
#   ggtitle(paste0("Spatial distribution of Seurat clusters at Res = 0.2 for ",
#                  gsub("_", " sample ", unique(seu@meta.data$sample_ident)[i]))
#           ) -> gpl[[i]]
# }
# 
# pdf(file = paste0(fp0, "spinal_cord_spatial_cluster_by_sample.pdf"))
# gpl
# dev.off()

################################# DE comparisons ###############################
library("dplyr")
library("readr")
library("writexl")
library("Seurat")
library("tibble")
library("ggplot2")
library("future")

## read saved prepared Seurat object
fp0 <- "/lustre1/project/stg_00104/vsc/Projects/Project_Theo/data/piol_data/RDS/"
seu <- readRDS(file = paste0(fp0, "SC_filtered_seurat.rds"))

## edit metadata table
seu@meta.data %>%
  dplyr::mutate(group = dplyr::case_when(
    orig.ident == "CtrlFUSH" ~ "Ctrl",
    orig.ident == "CtrlSOD1" ~ "Ctrl",
    orig.ident == "CtrlTDP43" ~ "Ctrl",
    orig.ident == "FUSH" ~ "FUSH",
    orig.ident == "SOD1" ~ "SOD1",
    orig.ident == "TDP43" ~ "TDP43")) %>%
  as.data.frame() -> seu@meta.data

## subset each experimental group into its own Seurat object
seu0 <- SetIdent(object = seu, value = "orig.ident")
sub1 <- subset(seu0, subset = orig.ident == c("CtrlFUSH"))
sub2 <- subset(seu0, subset = orig.ident == c("CtrlSOD1"))
sub3 <- subset(seu0, subset = orig.ident == c("CtrlTDP43"))
sub4 <- subset(seu0, subset = orig.ident == c("FUSH"))
sub5 <- subset(seu0, subset = orig.ident == c("SOD1"))
sub6 <- subset(seu0, subset = orig.ident == c("TDP43"))
merged1 <- merge(x = sub1, y = sub4) ## c("CtrlFUSH", "FUSH")
merged2 <- merge(x = sub3, y = sub6) ## c("CtrlTDP43", "TDP43")
merged3 <- merge(x = sub1, y = c(sub3, sub6)) ## c("CtrlFUSH", "CtrlTDP43", "TDP43")
merged4 <- merge(x = sub1, y = c(sub3, sub4)) ## c("CtrlFUSH", "CtrlTDP43", "FUSH")
merged5 <- merge(x = sub2, y = sub5) ## c("CtrlSOD1", "SOD1")
merged1 <- SetIdent(object = merged1, value = "orig.ident")
merged2 <- SetIdent(object = merged2, value = "orig.ident")
merged3 <- SetIdent(object = merged3, value = "group")
merged4 <- SetIdent(object = merged4, value = "group")
merged5 <- SetIdent(object = merged5, value = "orig.ident")
levels(merged1)
levels(merged2)
levels(merged3)
levels(merged4)
levels(merged5)

## Differential expression analysis
seu0 <- SetIdent(object = seu, value = "orig.ident")
seu0 %>%
  FindAllMarkers(test.use = "wilcox") %>%
  as.data.frame() %>%
  dplyr::arrange(p_val_adj) %>%
  dplyr::select(gene, cluster, everything()) %>%
  dplyr::rename(group = cluster) %>%
  as.data.frame() -> seu_mark0
rownames(seu_mark0) <- NULL
merged1 %>%
  FindAllMarkers(test.use = "wilcox") %>%
  as.data.frame() %>%
  dplyr::arrange(p_val_adj) %>%
  dplyr::select(gene, cluster, everything()) %>%
  dplyr::rename(group = cluster) %>%
  as.data.frame() -> seu_mark1
rownames(seu_mark1) <- NULL
# merged2 %>%
#   FindAllMarkers(test.use = "wilcox") %>%
#   as.data.frame() %>%
#   dplyr::arrange(p_val_adj) %>%
#   dplyr::select(gene, cluster, everything()) %>%
#   dplyr::rename(group = cluster) %>%
#   as.data.frame() -> seu_mark2
# rownames(seu_mark2) <- NULL
# merged3 %>%
#   FindAllMarkers(test.use = "wilcox") %>%
#   as.data.frame() %>%
#   dplyr::arrange(p_val_adj) %>%
#   dplyr::select(gene, cluster, everything()) %>%
#   dplyr::rename(group = cluster) %>%
#   as.data.frame() -> seu_mark3
# rownames(seu_mark3) <- NULL
merged4 %>%
  FindAllMarkers(test.use = "wilcox") %>%
  as.data.frame() %>%
  dplyr::arrange(p_val_adj) %>%
  dplyr::select(gene, cluster, everything()) %>%
  dplyr::rename(group = cluster) %>%
  as.data.frame() -> seu_mark4
rownames(seu_mark4) <- NULL
merged5 %>%
  FindAllMarkers(test.use = "wilcox") %>%
  as.data.frame() %>%
  dplyr::arrange(p_val_adj) %>%
  dplyr::select(gene, cluster, everything()) %>%
  dplyr::rename(group = cluster) %>%
  as.data.frame() -> seu_mark5
rownames(seu_mark5) <- NULL

## save results as Excel spreadsheets with each comparison on a different sheet
seu_markers <- list(seu_mark1,
                    # seu_mark2, seu_mark3,
                    seu_mark4,
                    seu_mark5, seu_mark0)
names(seu_markers) <- c("ctrlFUS vs FUS",
                        # "ctrlTDP vs TDP",
                        # "ctrlFUS and ctrlTDP vs TDP",
                        "ctrlFUS and ctrlTDP vs FUS",
                        "ctrlSOD vs SOD", "all vs all")
writexl::write_xlsx(path = paste0(fp0, "seurat_DE_res.xlsx"),
                    x = seu_markers)

#### obsolete ??

## read saved prepared Seurat object
# seu <- readRDS(file = paste0(fp2, "seurat_", pn, ".rds"))

## DE on clusters
# seu <- SetIdent(object = seu, value = "RNA_snn_res.0.1")
# levels(seu)
# seu %>%
#   FindAllMarkers(test.use = "wilcox") %>%
#   as.data.frame() %>% 
#   dplyr::arrange(p_val_adj) %>%
#   dplyr::select(gene, cluster, everything()) %>%
#   as.data.frame() -> seu_mark_c
# saveRDS(object = seu_mark_c, file = paste0(fp2, pn, "_res_0.1_de.rds"))