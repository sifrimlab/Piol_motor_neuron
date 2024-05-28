# For resolve:
# Analyze spinal cord DAPI as single cells and identify the motor neuron cluster
# You can try these markers: Chat/ Chodl, you can send me the umap/tsne with the
# expression pattern of these two markers and we will see if we can find the
# motor neuron cluster and start the differential analysis from there.

## load libraries
library("dplyr")
library("readr")
library("writexl")
library("Seurat")
library("tibble")
library("ggplot2")
library("future")
library("forcats")
library("gridExtra")

## project name
pn <- "spinal_cord"

## input file path
fp <- "~/Documents/tmp/dacruz/202305_piol/count_matrices_correct/SC_Resolve_Images/"

## out path
out <- "~/Documents/tmp/dacruz/202310_piol/data/RDS/"

## process Excel count matrices
# lf <- list.files(path = fp, pattern = "*.csv")
# all_dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
# all_dfs2 <- lapply(all_dfs1, function(i){i %>%
#     tibble::column_to_rownames(var = "cell_label")})
# all_dfs3 <- lapply(all_dfs2, function(i){t(i)})
# 
# ## create list of seurat objects from the transformed count matrices
# fls1 <- lapply(all_dfs3, function(i) CreateSeuratObject(counts = i,
#                                                         min.cells = 1,
#                                                         min.features = 1))
# fls2 <- list()
# for (i in 1:length(fls1)) {
#   ## get experimental control condition and add to metadata
#   fls1[[i]]@meta.data %>%
#     dplyr::mutate(orig.ident = gsub("_.*", "", lf)[i],
#                   sample_name = gsub("_labeled_count-matrix.csv", "", lf)[i]
#                   ) -> fls2[[i]]
#   fls1[[i]]@meta.data <- fls2[[i]]         ## add sample-wise metadata to object
#   fls1[[i]]@meta.data$sample <- paste(i)   ## add sample number
#   ## add x and y coords
#   fls1[[i]]@meta.data$y <- as.numeric(
#     gsub("row", "", gsub("_.*", "", rownames(fls2[[i]]))))
#   fls1[[i]]@meta.data$x <- as.numeric(
#     gsub("col", "", gsub(".*_", "", gsub("_[^_]+$", "", rownames(fls2[[i]])))))
# }
# 
# ## merge list of processed Excel count matrices into single Seurat object
# merge(x = fls1[[1]], y = fls1[2:length(fls1)]) %>%
#   NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>%
#   FindVariableFeatures(selection.method = "vst") %>%
#   ScaleData() %>%
#   RunPCA(npcs = 10, verbose = FALSE) %>%
#   FindNeighbors(reduction = "pca", dims = 1:10, compute.SNN = TRUE,
#                 nn.method = "rann", nn.eps = 0) %>%
#   FindClusters(resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)) %>%
#   RunUMAP(reduction = "pca", dims = 1:10) -> seu
# 
# ## save prepared Seurat object
# saveRDS(object = seu, file = paste0(out, "seurat_", pn, ".rds"))
# 
# ## make marker gene plots for GSE161621
# mn_genes <- c(
#   "Cdh6", #	fast-firing MN marker
#   "Chat", ## universal MN marker
#   "Chodl", # fast-firing MN marker
#   "Kcnq5", # fast-firing MN marker
#   "Kcnt2", # fast-firing MN marker
#   "Prkcb", # fast-firing MN marker
#   "Kcnd2", # slow-firing MN marker
#   "Prkcd", # slow-firing MN marker
#   "Sv2a"  # slow-firing MN marker
# )
## new marker gene list
# readxl::read_excel(path = paste0(out, "neural_marker_genes.xlsx")) %>%
#   dplyr::select(Gene) %>% 
#   dplyr::pull() -> mn_genes
# 
# saveRDS(object = mn_genes, file = paste0(out, "mn_genes.rds"))
# 
# read_tsv(file = paste0(out, "PanglaoDB_markers_27_Mar_2020.tsv")) %>% 
#   as.data.frame() -> df1
# 
# saveRDS(object = df1, file = paste0(out, "PanglaoDB_markers_27_Mar_2020.rds"))

## read RDS files
mn_genes <- readRDS(file = paste0(out, "mn_genes.rds"))
df1 <- readRDS(file = paste0(out, "PanglaoDB_markers_27_Mar_2020.rds"))
seu <- readRDS(file = paste0(out, "seurat_", pn, ".rds"))

## annotate cluster numbers with putative cell types
seu@meta.data %>% 
  dplyr::mutate(
    cluster_annotations = dplyr::case_when(
      RNA_snn_res.0.4 == 0 ~ "non-myelinating SCs?",
      RNA_snn_res.0.4 == 1 ~ "oligodendrocytes",
      RNA_snn_res.0.4 == 2 ~ "astrocytes?",
      RNA_snn_res.0.4 == 3 ~ "oligodendrocytes",
      RNA_snn_res.0.4 == 4 ~ "unknown",
      RNA_snn_res.0.4 == 5 ~ "OPC?",
      RNA_snn_res.0.4 == 6 ~ "unknown",
      RNA_snn_res.0.4 == 7 ~ "unknown",
      RNA_snn_res.0.4 == 8 ~ "microglia",
      RNA_snn_res.0.4 == 9 ~ "motor_neurons",
      RNA_snn_res.0.4 == 10 ~ "unknown",
      RNA_snn_res.0.4 == 11 ~ "astrocytes?")) %>%
  as.data.frame() -> seu@meta.data

## create PDF with cluster and expression UMAPs, dotplots and violin plots
# seu <- SetIdent(object = seu, value = "RNA_snn_res.0.4")
# seu <- SetIdent(object = seu, value = "cluster_annotations")
# seu <- SetIdent(object = seu, value = "RNA_snn_res.0.6")
# gg1 <- list(); gg2 <- list()
# for (i in 1:length(mn_genes)){
#   FeaturePlot(object = seu, features = mn_genes[[i]], reduction = "umap",
#               cols = c("azure3", "black"), pt.size = 1) -> gg1[[i]]
#   VlnPlot(object = seu, features = mn_genes[[i]], pt.size = 0) +
#     theme(axis.text.x = element_text(angle = 45, size = 10)) -> gg2[[i]]
# }
# pdf(file = paste0(out, "SC_spatial_marker_plots_annotated_v2_", Sys.Date(), ".pdf"))
# # pdf(file = paste0(out, "spinal_cord_new_marker_genes_v2.pdf"))
# UMAPPlot(seu, reduction = "umap", pt.size = 1, label = TRUE,
#          label.size = 5) + NoLegend() +
#   # labs(title = paste0("UMAP at Res 0.4")) +
#   labs(title = paste0("Tentative UMAP cluster annotations")) +
#   theme(plot.title = element_text(hjust = 0.5))
# DotPlot(seu, features = mn_genes) +
#   RotatedAxis() +
#   theme(axis.text.x = element_text(angle = 45, size = 8))
# for (i in 1:length(mn_genes)){
#   print(gg1[[i]])
#   print(gg2[[i]])
# }
# dev.off()

## spatial plot
# pdf(file = paste0(out, "spinal_cord_spatial_cluster_plots_", Sys.Date(), ".pdf"))
# gg3 <- list()
# for (i in 1:length(unique(seu@meta.data$sample_name))){
#   seu@meta.data %>%
#     dplyr::filter(sample_name == unique(seu@meta.data$sample_name)[i]) %>%
#     ggplot(aes(x = x, y = y, color = RNA_snn_res.0.4)) +
#     geom_point(size = 0.5) +
#     ggtitle(paste0(unique(seu@meta.data$sample_name)[i])) -> gg3[[i]]
# }
# print(gg3)
# dev.off()

## annotated spatial plots #1
seu <- SetIdent(object = seu, value = "cluster_annotations")
gg4 <- list()
for (i in 1:length(unique(seu@meta.data$sample_name))){
  seu@meta.data %>%
    dplyr::filter(sample_name == unique(seu@meta.data$sample_name)[i]) %>%
    ggplot(aes(x = x, y = y, color = cluster_annotations)) +
    geom_point(size = 0.5) +
    ggtitle(paste0(unique(seu@meta.data$sample_name)[i])) -> gg4[[i]]
}
# pdf(file = paste0(out, "SC_spatial_annotated_cluster_plots_", Sys.Date(), ".pdf"))
# print(gg4)
# dev.off()

## select annotated spatial plots #2 - SOD1
pdf(file = paste0(out, "SC_spatial_annotated_plots_SOD1_", Sys.Date(), ".pdf"))
gridExtra::grid.arrange(gg4[[3]] + theme(legend.position = "none"),
                        gg4[[11]] + theme(legend.position = "none"),
                        gg4[[4]] + theme(legend.position = "none"),
                        gg4[[12]] + theme(legend.position = "none"),
                        gg4[[5]] + theme(legend.position = "none"),
                        gg4[[13]] + theme(legend.position = "none"),
                        ncol = 2)
dev.off()

##############################################################################

# table(seu@meta.data$RNA_snn_res.0.4)

## read saved RDS object
# seu <- readRDS(file = paste0(out, "seurat_", pn, ".rds"))

## motor neuron datasets
# seu@meta.data %>%
#   dplyr::filter(RNA_snn_res.0.4 == 9) -> df1
# 
# as.data.frame(table(df1$sample_name)) %>%
#   dplyr::rename(sample_name = names(.)[1],
#                 cluster_9_obj_num = names(.)[2]) -> df2
# 
# as.data.frame(table(seu@meta.data$sample_name)) %>%
#   dplyr::rename(sample_name = names(.)[1],
#                 total_obj_num = names(.)[2]) -> df3
# 
# df2 %>% 
#   dplyr::left_join(df3, by = "sample_name") %>% 
#   dplyr::mutate(cluster_9_perc = cluster_9_obj_num / total_obj_num,
#                 group = as.factor(gsub("_.*", "", sample_name)),
#                 group = forcats::fct_relevel(
#                   group, c("CtrlFUSH", "CtrlTDP43", "TDP43", "FUSH",
#                            "CtrlSOD1", "SOD1"))) -> df4
# 
# df4 %>%
#   ggplot(aes(x = group, y = total_obj_num, color = group)) +
#   ylab("Total number of DAPI objects") +
#   geom_point() +
#   ggtitle("Total number of DAPI objects for each sample by group") -> b1
# 
# df4 %>%
#   ggplot(aes(x = group, y = cluster_9_obj_num, color = group)) +
#   ylab("Total number of DAPI objects in Cluster 9 (motor neurons)") +
#   geom_point() +
#   ggtitle("Total number of motor neuron DAPI objects for each sample by group") -> b2
# 
# df4 %>%
#   ggplot(aes(x = group, y = cluster_9_perc, color = group)) +
#   ylab("Percentage of motor neurons in total DAPI objects") +
#   geom_point() +
#   ggtitle(paste0("Percentage of motor neuron objects to total number of DAPI",
#                  " objects\nfor each sample by group")) -> b3
# 
# pdf(file = paste0(out, "spinal_cord_spatial_cluster_cell_percentage.pdf"))
# b1;b2;b3
# dev.off()

############# single cell differential expression - FindAllMarkers #############