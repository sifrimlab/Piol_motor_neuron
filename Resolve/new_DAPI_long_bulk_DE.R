## NEW RESOLVE BULK DE - DAPI long

## load libraries
library("dplyr")
library("readr")
library("tidyr")
library("tibble")
library("stringr")
library("ggplot2")
library("ggrepel")
library("writexl")
library("DESeq2")

## other parameters
options(scipen = 999)
fp0 <- "~/Documents/tmp/2023/202308_piol/new_plots/"

############################### SN1 "DAPI longitudinal" ########################
fp1 <- "~/Documents/tmp/2023/202305_piol/count_matrices_correct/SN_Resolve_MC_Images/DAPI/longitudinal/"
lf <- list.files(path = fp1, pattern = "*.csv")
all_dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp1, i))})
all_dfs2 <- lapply(all_dfs1, function(i){i %>%
    tibble::column_to_rownames(var = "cell_label")})
all_dfs3 <- lapply(all_dfs2, function(i){t(i)})
rl1 <- list()
for (i in 1:length(all_dfs3)) {
  rl1[[i]] <- as.data.frame(rowSums(all_dfs3[[i]]))
}

tibble::rownames_to_column(rl1[[1]], var = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[2]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[3]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[4]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[5]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[6]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[7]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[8]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[9]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[10]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[11]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[12]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[13]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[14]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[15]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[16]], var = "genes"), by = "genes") %>%
  as.data.frame() -> c1
c1[is.na(c1)] <- 0 ## coerce NA values to zero
names(c1) <- c("genes", paste0(gsub("-", "_", gsub("_labeled_count-matrix.csv", "", lf))))

## make metadata table
data.frame(file_name = as.factor(lf),
           sample = paste0(gsub("-", "_", gsub("_labeled_count-matrix.csv", "", lf))),
           stain = "DAPI", 
           slice = "transversal",
           stain_slice = "DAPI_transversal",
           batch = "SN_Resolve_1",
           group = as.factor(gsub("_.*", "", lf)),
           control = as.factor(c(rep("CtrlFUSH+CtrlTDP43", 2), rep("CtrlSOD1", 3),
                                 rep("CtrlFUSH+CtrlTDP43", 2), rep("FUSH", 3),
                                 rep("SOD1", 3), rep("TDP43", 3)))) -> m1
rownames(m1) <- m1$sample

############################### SN2 "DAPI transversal" #########################
fp1 <- "~/Documents/tmp/2023/202305_piol/count_matrices_correct/SN_Resolve_Images_2/DAPI/longitudinal/"
lf <- list.files(path = fp1, pattern = "*.csv")
all_dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp1, i))})
all_dfs2 <- lapply(all_dfs1, function(i){i %>%
    tibble::column_to_rownames(var = "cell_label")})
all_dfs3 <- lapply(all_dfs2, function(i){t(i)})
rl1 <- list()
for (i in 1:length(all_dfs3)) {
  rl1[[i]] <- as.data.frame(rowSums(all_dfs3[[i]]))
}

tibble::rownames_to_column(rl1[[1]], var = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[2]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[3]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[4]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[5]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[6]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[7]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[8]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[9]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[10]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[11]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[12]], var = "genes"), by = "genes") %>%
  dplyr::full_join(tibble::rownames_to_column(rl1[[13]], var = "genes"), by = "genes") %>%
  as.data.frame() -> c2
c2[is.na(c2)] <- 0 ## coerce NA values to zero
names(c2) <- c("genes", paste0(gsub("-", "_", gsub("_labeled_count-matrix.csv", "", lf))))

## make metadata table
data.frame(file_name = as.factor(lf),
           sample = paste0(gsub("-", "_", gsub("_labeled_count-matrix.csv", "", lf))),
           stain = "DAPI", 
           slice = "transversal",
           stain_slice = "DAPI_transversal",
           batch = "SN_Resolve_2",
           group = as.factor(gsub("_.*", "", lf)),
           control = as.factor(c(rep("CtrlFUSH+CtrlTDP43", 2), rep("CtrlSOD1", 3),
                                 "CtrlFUSH+CtrlTDP43", rep("FUSH", 3),
                                 rep("SOD1", 2), rep("TDP43", 2)))
           ) -> m2
rownames(m2) <- m2$sample

## combine SN1 and SN2
c1 %>% 
  dplyr::full_join(c2, by = "genes") %>%
  dplyr::filter(!grepl("FP", genes)) %>%
  tibble::column_to_rownames(var = "genes") %>% 
  as.data.frame() -> SN_dapi_all
SN_dapi_all[is.na(SN_dapi_all)] <- 0 ## coerce NA values to zero

dplyr::bind_rows(m1, m2) %>% 
  as.data.frame() -> sn_meta_dapi

############################### SN2 "DAPI transversal" #########################

## run DESeq2 on "group"
DESeqDataSetFromMatrix(countData = as.matrix(SN_dapi_all),
                       colData = sn_meta_dapi,
                       design = ~ control) -> dds1
keep1 <- rowSums(counts(dds1)) > 1
dds1 <- DESeq(dds1[keep1, ])
# resultsNames(dds1)
rld1 <- rlogTransformation(dds1)
PCA1 <- plotPCA(rld1, intgroup = "control") +
        theme(plot.title = element_text(hjust = 0.5)) +
        # geom_text_repel(aes(label = colnames(rld1)),
        #                 arrow = arrow(length = unit(0.03, "npc"),
        #                 type = "closed", ends = "first"), force = 5) +
        ggtitle(paste0("Principal component analysis of SN1 and SN2 pseudobulk counts"))
PCA2 <- plotPCA(rld1, intgroup = "batch") +
        theme(plot.title = element_text(hjust = 0.5)) +
        # geom_text_repel(aes(label = colnames(rld1)),
        #                 arrow = arrow(length = unit(0.03, "npc"),
        #                 type = "closed", ends = "first"), force = 5) +
        ggtitle(paste0("Principal component analysis of SN1 and SN2 pseudobulk counts"))

## DE comparisons control FUS vs FUS
res1 <- results(dds1, contrast = c("control", "CtrlFUSH+CtrlTDP43", "FUSH"), alpha = 0.05)
merge(as.data.frame(res1), as.data.frame(counts(dds1, normalized = TRUE)),
      by = "row.names", sort = FALSE) %>%
  as.data.frame() %>%
  dplyr::rename(genes = names(.)[1]) %>%
  dplyr::arrange(padj) -> rd1

## DE comparisons control TDP vs TDP
res2 <- results(dds1, contrast = c("control", "CtrlFUSH+CtrlTDP43", "TDP43"), alpha = 0.05)
merge(as.data.frame(res2), as.data.frame(counts(dds1, normalized = TRUE)),
      by = "row.names", sort = FALSE) %>%
  as.data.frame() %>%
  dplyr::rename(genes = names(.)[1]) %>%
  dplyr::arrange(padj) -> rd2

## DE comparisons control SOD vs SOD
res3 <- results(dds1, contrast = c("control", "CtrlSOD1", "SOD1"), alpha = 0.05)
merge(as.data.frame(res3), as.data.frame(counts(dds1, normalized = TRUE)),
      by = "row.names", sort = FALSE) %>%
  as.data.frame() %>%
  dplyr::rename(genes = names(.)[1]) %>%
  dplyr::arrange(padj) -> rd3

rd3 %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::mutate(threshold = as.factor(abs(log2FoldChange) > 1 & padj < 0.05)) %>%
  dplyr::select(genes:padj, threshold) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
  theme(legend.position = "none") +
  geom_point(alpha = 0.8, size = 0.8) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_vline(xintercept = -1, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
  xlab("log2 fold change") + ylab("-log10 padj-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text_repel(aes(label = genes)) +
  scale_color_manual(values = c("black", "red")) +
  ggtitle(paste0("DE genes between CtrlSOD1 and SOD1 in DAPI longitudinal")) -> v1


## save pseudobulk_gene_counts
pdf(file = paste0(fp0, "SN1_SN2_DAPI_longitudinal_pseudobulk_DE_volcano.pdf"))
print(v1)
dev.off()

## save results
SN_dapi_all %>%
  tibble::rownames_to_column(var = "genes") %>%
  dplyr::arrange(genes) %>%
  base::as.data.frame() -> y1
res_list <- list(y1, rd1, rd2, rd3)
names(res_list) <- c("raw pseudobulk segmented counts",
                     "ctrlFUSH and ctrlTDP43 vs FUSH",
                     "ctrlFUSH and ctrlTDP43 vs TDP43",
                     "ctrlSOD1 vs SOD1")
writexl::write_xlsx(x = res_list, path = paste0(
  fp0, "SN1_SN2_DAPI_longitudinal_pseudobulk_DESeq2_res.xlsx"))

## save pseudobulk_gene_counts
pdf(file = paste0(fp0, "SN1_SN2_DAPI_longitudinal_pseudobulk_PCAs.pdf"))
print(PCA1)
print(PCA2)
dev.off()
