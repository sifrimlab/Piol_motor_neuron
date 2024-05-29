## NEW SEGMENTATION VS. MASK. We are comparing the segmented objects in
## NFH to everything else in the tissue (i.e. everything that is outside of NFH
## and DAPI), but only in the only Control TDP43 and Control Fus samples. This
## script will collate the counts from the conditions mentioned above into a
## pseudobulk matrix and perform DE using DESeq2 and then we need to filter the
## collated count matrix to contain only Control TDP43 and Control Fus samples!!

## load libraries
library("dplyr")
library("readr")
library("writexl")
library("tidyr")
library("tibble")
library("stringr")
library("forcats")
library("ggplot2")
library("ggrepel")
library("DESeq2")

## other parameters
options(scipen = 999)
fp0 <- "~/Documents/tmp/2023/202308_piol/data/"

##################################### SN1 MASK #################################
fp <- paste0("~/Documents/tmp/2023/202308_piol/data/new_segmentation_vs_mask/",
             "SN_Resolve_MC_Images/")
lf <- list.files(path = fp, pattern = "*.csv")
dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    tibble::column_to_rownames(var = "cell_label")})
dfs3 <- lapply(dfs2, function(i){t(i)})
rl <- list()
for (i in 1:length(dfs3)) {
  rl[[i]] <- as.data.frame(rowSums(dfs3[[i]]))
}
tibble::rownames_to_column(rl[[1]], var = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[2]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[3]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[4]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[5]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[6]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[7]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[8]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[9]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[10]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[11]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[12]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[13]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[14]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[15]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[16]], var = "genes"), by = "genes") %>%
  dplyr::filter(!grepl("FP", genes)) %>% ## throw out false positives
  dplyr::arrange(desc(genes)) %>%
  # tibble::column_to_rownames(var = "genes") %>%
  as.data.frame() -> cn1
cn1[is.na(cn1)] <- 0 ## coerce NA values to zero
names(cn1) <- c("genes", paste0(
  gsub("-", "_", gsub("_DAPI_labeled_noDAPI-noNFH-negative-count-matrix.csv", "", lf)
       ), "_SN1_mask"))

##################################### SN2 MASK #################################
fp <- paste0("~/Documents/tmp/2023/202308_piol/data/new_segmentation_vs_mask/",
             "SN_Resolve_Images_2/")
lf <- list.files(path = fp, pattern = "*.csv")
dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    tibble::column_to_rownames(var = "cell_label")})
dfs3 <- lapply(dfs2, function(i){t(i)})
rl <- list()
for (i in 1:length(dfs3)) {
  rl[[i]] <- as.data.frame(rowSums(dfs3[[i]]))
}
tibble::rownames_to_column(rl[[1]], var = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[2]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[3]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[4]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[5]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[6]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[7]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[8]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[9]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[10]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[11]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[12]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[13]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[14]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[15]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[16]], var = "genes"), by = "genes") %>%
  dplyr::filter(!grepl("FP", genes)) %>% ## throw out false positives
  dplyr::arrange(desc(genes)) %>%
  # tibble::column_to_rownames(var = "genes") %>%
  as.data.frame() -> cn2
cn2[is.na(cn2)] <- 0 ## coerce NA values to zero
names(cn2) <- c("genes", paste0(
  gsub("-", "_", gsub("_DAPI_labeled_noDAPI-noNFH-negative-count-matrix.csv", "", lf)
       ), "_SN2_mask"))

############################# "SN1 NFH_transversal" ############################
fp <- paste0("~/Documents/tmp/2023/202305_piol/count_matrices_correct/with_area/",
             "SN_Resolve_MC_Images/non_expanded/NFH/")
lf <- list.files(path = fp, pattern = "*.csv")
dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    dplyr::filter(area < 2500) %>%
    tibble::column_to_rownames(var = "cell_label")})
dfs3 <- lapply(dfs2, function(i){t(i)})
rl <- list()
for (i in 1:length(dfs3)) {
  rl[[i]] <- as.data.frame(rowSums(dfs3[[i]]))
}
tibble::rownames_to_column(rl[[1]], var = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[2]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[3]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[4]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[5]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[6]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[7]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[8]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[9]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[10]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[11]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[12]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[13]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[14]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[15]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[16]], var = "genes"), by = "genes") %>%
  dplyr::filter(!grepl("FP", genes)) %>% ## throw out false positives
  dplyr::arrange(desc(genes)) %>%
  # tibble::column_to_rownames(var = "genes") %>%
  as.data.frame() -> cp1
cp1[is.na(cp1)] <- 0 ## coerce NA values to zero
names(cp1) <- c("genes", paste0(
  gsub("-", "_", gsub("_labeled_non-expanded_count-matrix_with-area.csv", "", lf)), "_SN1_seg"))

############################# "SN2 NFH_transversal" ############################
fp <- paste0("~/Documents/tmp/2023/202305_piol/count_matrices_correct/with_area/",
             "SN_Resolve_Images_2/non_expanded/NFH/")
lf <- list.files(path = fp, pattern = "*.csv")
dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    dplyr::filter(area < 2500) %>%
    tibble::column_to_rownames(var = "cell_label")})
dfs3 <- lapply(dfs2, function(i){t(i)})
rl <- list()
for (i in 1:length(dfs3)) {
  rl[[i]] <- as.data.frame(rowSums(dfs3[[i]]))
}
tibble::rownames_to_column(rl[[1]], var = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[2]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[3]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[4]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[5]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[6]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[7]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[8]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[9]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[10]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[11]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[12]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[13]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[14]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[15]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl[[16]], var = "genes"), by = "genes") %>%
  dplyr::filter(!grepl("FP", genes)) %>% ## throw out false positives
  dplyr::arrange(desc(genes)) %>%
  # tibble::column_to_rownames(var = "genes") %>%
  as.data.frame() -> cp2
cp2[is.na(cp2)] <- 0 ## coerce NA values to zero
names(cp2) <- c("genes", paste0(
  gsub("-", "_", gsub("_labeled_non-expanded_count-matrix_with-area.csv", "", lf)), "_SN2_seg"))

################################# COLLATE THE MATRICES #########################
cp1 %>%
  dplyr::mutate(total_counts = rowSums(across(where(is.numeric)))) %>%
  dplyr::select(genes, total_counts, everything()) %>%
  dplyr::arrange(genes) %>%
  as.data.frame() -> cp1s
cp2 %>%
  dplyr::mutate(total_counts = rowSums(across(where(is.numeric)))) %>%
  dplyr::select(genes, total_counts, everything()) %>%
  dplyr::arrange(genes) %>%
  as.data.frame() -> cp2s

cn1 %>%
  dplyr::mutate(total_counts = rowSums(across(where(is.numeric)))) %>%
  dplyr::select(genes, total_counts, everything()) %>%
  dplyr::arrange(genes) %>%
  as.data.frame() -> cn1s
cn2 %>%
  dplyr::mutate(total_counts = rowSums(across(where(is.numeric)))) %>%
  dplyr::select(genes, total_counts, everything()) %>%
  dplyr::arrange(genes) %>%
  as.data.frame() -> cn2s

cp1s %>%
  dplyr::select(-total_counts) %>%
  dplyr::filter(genes != "area") %>%
  dplyr::full_join((cp2s %>%
                       dplyr::select(-total_counts) %>%
                       dplyr::filter(genes != "area")) , by = "genes") %>%
    dplyr::full_join((cn1s %>%
                       dplyr::select(-total_counts)) , by = "genes") %>%
    dplyr::full_join((cn2s %>%
                       dplyr::select(-total_counts)) , by = "genes") %>%
  dplyr::select(genes, contains("CtrlFUSH"), contains("CtrlTDP43")) %>%
  tibble::column_to_rownames(var = "genes") %>%
  as.data.frame() -> c_df
c_df[is.na(c_df)] <- 0 ## coerce NA values to zero

## build excel tables
# cl1 <- list(cp1s, cp2s, cn1s, cn2s)
# names(cl1) <- c("SN1_NFH_segmented_objects",
#                 "SN2_NFH_segmented_objects",
#                 "SN1_NFH_mask",
#                 "SN2_NFH_mask")
# writexl::write_xlsx(x = cl1, path = paste0(
#   fp0, "SN_NFH_segmented_objects_and_mask_pseudobulk_cts.xlsx"))

## save count matrices
# saveRDS(object = c_df, file = paste0(fp0, "SN_segmented_objects_and_mask_counts.rds"))
# c_df <- readRDS(file = paste0(fp0, "SN_segmented_objects_and_mask_counts.rds"))
# 
c_df %>%
  t() %>% as.data.frame() %>%
  tibble::rownames_to_column(var = "sample") %>%
  dplyr::select(sample) %>%
    dplyr::mutate(group = dplyr::case_when(
    stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH",
    # stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
    stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlTDP43" #,
    # stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
    # stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
    # stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43"
    )) %>%
  dplyr::mutate(control = dplyr::case_when(
    stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
    # stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
    stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43" #,
    # stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
    # stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
    # stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43"
    )) %>%
  dplyr::mutate(batch = dplyr::case_when(
    stringr::str_detect(sample, pattern = "SN1") ~ "SN1",
    stringr::str_detect(sample, pattern = "SN2") ~ "SN2")) %>%
  dplyr::mutate(segmentation = dplyr::case_when(
    stringr::str_detect(sample, pattern = "seg") ~ "segmented",
    stringr::str_detect(sample, pattern = "mask") ~ "mask")) %>%
  # dplyr::mutate(control = as.factor(control),
  #               control = forcats::fct_relevel(control,
  #                                              c("CtrlFUSH+CtrlTDP43",
  #                                                "FUSH", "TDP43",
  #                                                "CtrlSOD1", "SOD1"))) %>%
  dplyr::mutate_all(as.factor) %>%
  as.data.frame() -> meta_data
rownames(meta_data) <- meta_data$sample

## run DESeq2 on group + batch + segmentation
DESeqDataSetFromMatrix(countData = as.matrix(c_df),
                       colData = meta_data,
                       design = ~ group + batch + segmentation) -> dds1
keep1 <- rowSums(counts(dds1)) > 1
dds1 <- DESeq(dds1[keep1, ])
# 
# ## run DESeq2 on control + batch + segmentation
# DESeqDataSetFromMatrix(countData = as.matrix(c_df),
#                        colData = meta_data,
#                        design = ~ control + batch + segmentation) -> dds2
# keep2 <- rowSums(counts(dds2)) > 1
# dds2 <- DESeq(dds2[keep2, ])
# 
# # resultsNames(dds)
rld1 <- rlogTransformation(dds1)
PCA1 <- plotPCA(rld1, intgroup = "segmentation") +
  theme(plot.title = element_text(hjust = 0.5)) +
  # geom_text_repel(aes(label = colnames(rld1)),
  #                 arrow = arrow(length = unit(0.03, "npc"),
  #                 type = "closed", ends = "first"), force = 5) +
  ggtitle(paste0("PCA of SN NFH segmented objects and mask bulk counts"))
PCA2 <- plotPCA(rld1, intgroup = "batch") +
  theme(plot.title = element_text(hjust = 0.5)) +
  # geom_text_repel(aes(label = colnames(rld1)),
  #                 arrow = arrow(length = unit(0.03, "npc"),
  #                 type = "closed", ends = "first"), force = 5) +
  ggtitle(paste0("PCA of SN NFH segmented objects and mask bulk counts"))
# PCA3 <- plotPCA(rld1, intgroup = "control") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   # geom_text_repel(aes(label = colnames(rld)),
#   #                 arrow = arrow(length = unit(0.03, "npc"),
#   #                 type = "closed", ends = "first"), force = 5) +
#   ggtitle(paste0("PCA of SN NFH segmented objects and mask bulk counts (control)"))

# ## size factor plots
# barplot(dds1$sizeFactor, main = "Size Factors for samples", col = "lavender",
#         cex.axis = 0.7, las = 2, ylab = "Size Factors", cex.names = 0.3) -> b1
 
# ############################### DE comparisons #################################

## DE comparisons control segmentation vs mask
res0 <- results(dds1, contrast = c("segmentation", "segmented", "mask"), alpha = 0.05)
merge(as.data.frame(res0), as.data.frame(counts(dds1, normalized = TRUE)),
      by = "row.names", sort = FALSE) %>%
  as.data.frame() %>%
  dplyr::rename(genes = names(.)[1]) %>%
  dplyr::arrange(padj) -> rd0

## DE comparisons control FUS vs FUS
# res1 <- results(dds1, contrast = c("group", "CtrlFUSH", "FUSH"), alpha = 0.05)
# merge(as.data.frame(res1), as.data.frame(counts(dds1, normalized = TRUE)),
#       by = "row.names", sort = FALSE) %>%
#   as.data.frame() %>%
#   dplyr::rename(genes = names(.)[1]) %>%
#   dplyr::arrange(padj) -> rd1
#
# ## DE comparisons control TDP vs TDP
# res2 <- results(dds1, contrast = c("group", "CtrlTDP43", "TDP43"), alpha = 0.05)
# merge(as.data.frame(res2), as.data.frame(counts(dds1, normalized = TRUE)),
#       by = "row.names", sort = FALSE) %>%
#   as.data.frame() %>%
#   dplyr::rename(genes = names(.)[1]) %>%
#   dplyr::arrange(padj) -> rd2
#
# ## DE comparisons control SOD vs SOD
# res3 <- results(dds1, contrast = c("group", "CtrlSOD1", "SOD1"), alpha = 0.05)
# merge(as.data.frame(res3), as.data.frame(counts(dds1, normalized = TRUE)),
#       by = "row.names", sort = FALSE) %>%
#   as.data.frame() %>%
#   dplyr::rename(genes = names(.)[1]) %>%
#   dplyr::arrange(padj) -> rd3
#
# ## DE comparisons control FUS & control TDP vs TDP
# res4 <- results(dds2, contrast = c("control", "TDP43", "CtrlFUSH+CtrlTDP43"), alpha = 0.05)
# merge(as.data.frame(res4), as.data.frame(counts(dds2, normalized = TRUE)),
#       by = "row.names", sort = FALSE) %>%
#   as.data.frame() %>%
#   dplyr::rename(genes = names(.)[1]) %>%
#   dplyr::arrange(padj) -> rd4
#
# ## DE comparisons control FUS & control TDP vs FUS
# res5 <- results(dds2, contrast = c("control", "FUSH", "CtrlFUSH+CtrlTDP43"), alpha = 0.05)
# merge(as.data.frame(res5), as.data.frame(counts(dds2, normalized = TRUE)),
#       by = "row.names", sort = FALSE) %>%
#   as.data.frame() %>%
#   dplyr::rename(genes = names(.)[1]) %>%
#   dplyr::arrange(padj) -> rd5

# ## save results and wide count matrix
c_df %>%  tibble::rownames_to_column(var = "genes") -> y1
res_list <- list(y1, rd0
                 #, rd1, rd2, rd3, rd4, rd5
                 )
names(res_list) <- c("raw pseudobulk counts",
                     "new_segmentation_vs_mask_DE"#,
                     # "ctrlFUSH vs FUSH",
                     # "ctrlTDP43 vs TDP43", "ctrlSOD1 vs SOD1",
                     # "ctrlFUSH and ctrlTDP43 vs TDP43",
                     # "ctrlFUSH and ctrlTDP43 vs FUSH"
                     )

writexl::write_xlsx(x = res_list, path = paste0(
  fp0, "new_SN_NFH_segmented_objects_vs_mask_pseudobulk_DE.xlsx"))

rd0 %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::mutate(threshold = as.factor(abs(log2FoldChange) > 1 & padj < 0.05)) %>%
  dplyr::select(genes:padj, threshold) %>%
  as.data.frame() -> res_vol

res_vol %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  geom_text_repel(aes(label = genes)) +
  theme(legend.position = "none") +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_vline(xintercept = -1, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
  xlab("log2 fold change") +
  ylab("-log10 padj-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("black", "red")) +
  ggtitle("NFH segmented objects vs mask (not NFH or DAPI)") -> v1

rd0 %>%
  dplyr::slice(1:5) %>%
  dplyr::select(genes, tidyr::contains("Ctrl")) %>%
  tidyr::gather(key = "sample", value = "count", -genes) %>%
  dplyr::mutate(segmentation = dplyr::case_when(
    stringr::str_detect(sample, pattern = "seg") ~ "segmented",
    stringr::str_detect(sample, pattern = "mask") ~ "mask"),
                genes = factor(genes, levels = unique(genes))) %>%
  dplyr::arrange(segmentation) %>%
  dplyr::mutate(segmentation = factor(segmentation, levels = unique(segmentation)),
                sample = factor(sample, levels = unique(sample))) %>%
  as.data.frame() %>%
  ggplot(aes(x = sample, y = log10(count), color = segmentation)) +
  geom_point() +
  theme(axis.text.x = element_blank()) +
  ggtitle(paste0("Gene counts of top 5 DE genes")) +
  facet_grid(~genes) -> t1

## save plots
pdf(file = paste0(fp0, "new_SN_NFH_segmented_objects_and_mask_pseudobulk_plots.pdf"))
PCA1; PCA2; v1; t1
dev.off()