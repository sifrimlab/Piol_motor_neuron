# The bar graphs should show the control TDP-43 + control FUS and FUS samples,
# and separate ones for control TDP-43 + control FUS and TDP-43 samples. This is
# an example of the bar graphs I would like to see:

# The datasets I would like to see these bar graphs are pseudobulk analyses from:
# NFH transversal
# DAPI transversal
# DAPI longitudinal
# ChAT+ transversal (NFH area > 1500)
# Mask of NFH+ vs everything that is not DAPI

## load libraries
library("dplyr")
library("ggplot2")
library("readr")
library("writexl")
library("tibble")
library("stringr")
library("forcats")

## plot output
fp0 <- "~/Documents/tmp/2023/202308_piol/new_plots/"

################################################################################
################ NFH 1500 gene count plots - "Chat pos area plots" #############
################################################################################

######################### "SN1 NFH_trans" CHAT pos (> 1500) ####################
fp <- paste0("~/Documents/tmp/2023/202305_piol/count_matrices_correct/with_area/",
             "SN_Resolve_MC_Images/non_expanded/NFH/")
lf <- list.files(path = fp, pattern = "*.csv")
dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    tibble::column_to_rownames(var = "cell_label") %>% 
    dplyr::filter(area > 1500) %>% 
    dplyr::select(-c(area, contains("FP")))})
dfs3 <- lapply(dfs2, function(i){t(i)})
rl1 <- list()
for (i in 1:length(dfs3)) {
  rl1[[i]] <- as.data.frame(rowSums(dfs3[[i]]))
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
  as.data.frame() -> cp1
cp1[is.na(cp1)] <- 0 ## coerce NA values to zero
names(cp1) <- c("genes", paste0(gsub("-", "_", gsub(
  "_labeled_non-expanded_count-matrix_with-area.csv", "", lf)), "_SN1_Chat_pos"))

######################### "SN2 NFH_trans" CHAT pos (> 1500) ####################
fp <- paste0("~/Documents/tmp/2023/202305_piol/count_matrices_correct/with_area/",
             "SN_Resolve_Images_2/non_expanded/NFH/")
lf <- list.files(path = fp, pattern = "*.csv")
dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    tibble::column_to_rownames(var = "cell_label") %>% 
    dplyr::filter(area > 1500) %>% 
    dplyr::select(-c(area, contains("FP")))})
dfs3 <- lapply(dfs2, function(i){t(i)})
rl1 <- list()
for (i in 1:length(dfs3)) {
  rl1[[i]] <- as.data.frame(rowSums(dfs3[[i]]))
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
  as.data.frame() -> cp2
cp2[is.na(cp2)] <- 0 ## coerce NA values to zero
names(cp2) <- c("genes", paste0(gsub("-", "_", gsub(
  "_labeled_non-expanded_count-matrix_with-area.csv", "", lf)), "_SN2_Chat_pos"))

######################### "SN1 NFH_trans" CHAT neg (< 1500) ####################
fp <- paste0("~/Documents/tmp/2023/202305_piol/count_matrices_correct/with_area/",
             "SN_Resolve_MC_Images/non_expanded/NFH/")
lf <- list.files(path = fp, pattern = "*.csv")
dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    tibble::column_to_rownames(var = "cell_label") %>% 
    dplyr::filter(area < 1500) %>% 
    dplyr::select(-c(area, contains("FP")))})
dfs3 <- lapply(dfs2, function(i){t(i)})
rl1 <- list()
for (i in 1:length(dfs3)) {
  rl1[[i]] <- as.data.frame(rowSums(dfs3[[i]]))
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
  as.data.frame() -> cn1
cn1[is.na(cn1)] <- 0 ## coerce NA values to zero
names(cn1) <- c("genes", paste0(gsub("-", "_", gsub(
  "_labeled_non-expanded_count-matrix_with-area.csv", "", lf)), "_SN1_Chat_neg"))

######################### "SN2 NFH_trans" CHAT neg (< 1500) ####################
fp <- paste0("~/Documents/tmp/2023/202305_piol/count_matrices_correct/with_area/",
             "SN_Resolve_Images_2/non_expanded/NFH/")
lf <- list.files(path = fp, pattern = "*.csv")
dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    tibble::column_to_rownames(var = "cell_label") %>% 
    dplyr::filter(area < 1500) %>% 
    dplyr::select(-c(area, contains("FP")))})
dfs3 <- lapply(dfs2, function(i){t(i)})
rl1 <- list()
for (i in 1:length(dfs3)) {
  rl1[[i]] <- as.data.frame(rowSums(dfs3[[i]]))
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
  as.data.frame() -> cn2
cn2[is.na(cn2)] <- 0 ## coerce NA values to zero
names(cn2) <- c("genes", paste0(gsub("-", "_", gsub(
  "_labeled_non-expanded_count-matrix_with-area.csv", "", lf)), "_SN2_Chat_neg"))

cp1 %>% 
  dplyr::left_join(cp2, by = "genes") %>%
  dplyr::left_join(cn1, by = "genes") %>% 
  dplyr::left_join(cn2, by = "genes") %>%
  dplyr::arrange(genes) -> NFH_1500_total_counts

NFH_1500_total_counts %>% 
  tidyr::gather(key = "sample", value = "counts", -genes) %>% 
  dplyr::mutate(
    motor_neuron = dplyr::case_when(
      grepl("Chat_pos", sample) ~ "Chat_pos",
      grepl("Chat_neg", sample) ~ "Chat_neg"),
    group = dplyr::case_when(
      stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH",
      stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
      stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlTDP43",
      stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
      stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
      stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43"),
    control = dplyr::case_when(
      stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
      stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
      stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
      stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
      stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
      stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43"),
    batch = dplyr::case_when(
      stringr::str_detect(sample, pattern = "SN1") ~ "SN1",
      stringr::str_detect(sample, pattern = "SN2") ~ "SN2"),
    control = as.factor(control),
    control = forcats::fct_relevel(control, c("CtrlFUSH+CtrlTDP43", "FUSH",
                                              "TDP43", "CtrlSOD1", "SOD1")),
    motor_neuron = forcats::fct_relevel(motor_neuron, c("CtrlFUSH+CtrlTDP43", "FUSH",
                                              "TDP43", "CtrlSOD1", "SOD1"))
    )  -> NFH_1500_total_counts_long

fgl1 <- list()
for (i in 1:length(unique(NFH_1500_total_counts_long$genes))) {
  NFH_1500_total_counts_long %>%
    dplyr::filter(genes == sort(unique(NFH_1500_total_counts_long$genes))[i]) %>%
    ggplot(aes(x = control, y = log(counts), fill = control)) +
    geom_boxplot() +
    geom_point(size = 2, color = "black") +
    ggtitle(paste0(unique(NFH_1500_total_counts_long$genes)[i],
                   " non-normalized pseudobulk counts")) +
    theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1),
          legend.position = "none") +
    facet_grid(~motor_neuron) -> fgl1[[i]]
}

## save pseudobulk_gene_counts
pdf(file = paste0(fp0, "SN_NFH_1500_Chat_pos_neg_pseudobulk_gene_count_boxplots.pdf"))
print(fgl1)
dev.off()

################################################################################
################ New Segmentation vs Mask - NFH vs "not NFH, not DAPI" #########
################################################################################

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
  as.data.frame() -> cp1
cp1[is.na(cp1)] <- 0 ## coerce NA values to zero
names(cp1) <- c("genes", paste0(
  gsub("-", "_", gsub("_labeled_non-expanded_count-matrix_with-area.csv", "", lf)),
  "_SN1_seg"))

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
  as.data.frame() -> cp2
cp2[is.na(cp2)] <- 0 ## coerce NA values to zero
names(cp2) <- c("genes", paste0(
  gsub("-", "_", gsub("_labeled_non-expanded_count-matrix_with-area.csv", "", lf)),
  "_SN2_seg"))

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
  # dplyr::select(genes, contains("CtrlFUSH"), contains("CtrlTDP43")) %>%
  as.data.frame() -> new_seg_df
new_seg_df[is.na(new_seg_df)] <- 0 ## coerce NA values to zero

new_seg_df %>% 
  tidyr::gather(key = "sample", value = "counts", -genes) %>% 
  dplyr::mutate(
    segmentation = dplyr::case_when(
      stringr::str_detect(sample, pattern = "seg") ~ "seg",
      stringr::str_detect(sample, pattern = "mask") ~ "mask"),
    group = dplyr::case_when(
      stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH",
      stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
      stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlTDP43",
      stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
      stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
      stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43"),
    control = dplyr::case_when(
      stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
      stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
      stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
      stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
      stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
      stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43"),
    batch = dplyr::case_when(
      stringr::str_detect(sample, pattern = "SN1") ~ "SN1",
      stringr::str_detect(sample, pattern = "SN2") ~ "SN2"),
    control = as.factor(control),
    control = forcats::fct_relevel(control, c("CtrlFUSH+CtrlTDP43", "FUSH",
                                              "TDP43", "CtrlSOD1", "SOD1"))) %>% 
  dplyr::arrange(genes) -> new_seg_df_long

fgl2 <- list()
for (i in 1:length(unique(new_seg_df_long$genes))) {
  new_seg_df_long %>%
    dplyr::filter(genes == sort(unique(new_seg_df_long$genes))[i]) %>%
    ggplot(aes(x = control, y = log(counts), fill = control)) +
    geom_boxplot() +
    geom_point(size = 2, color = "black") +
    ggtitle(paste0(unique(new_seg_df_long$genes)[i],
                   " non-normalized pseudobulk counts")) +
    theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1),
          legend.position = "none") +
    facet_grid(~segmentation) -> fgl2[[i]]
}

## save pseudobulk_gene_counts
pdf(file = paste0(fp0, "SN_seg_vs_mask_pseudobulk_gene_count_boxplots.pdf"))
print(fgl2)
dev.off()

################################################################################
########################### DAPI transversal SN ################################
################################################################################

############################### SN1 "DAPI transversal" #########################
fp1 <- "~/Documents/tmp/2023/202305_piol/count_matrices_correct/SN_Resolve_MC_Images/DAPI/transversal/"
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
names(c1) <- c("genes", paste0(gsub("-", "_", gsub("_labeled_count-matrix.csv", "", lf)),
                               "_SN1"))

############################### SN2 "DAPI transversal" #########################
fp1 <- "~/Documents/tmp/2023/202305_piol/count_matrices_correct/SN_Resolve_Images_2/DAPI/transversal/"
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
  as.data.frame() -> c2
c2[is.na(c2)] <- 0 ## coerce NA values to zero
names(c2) <- c("genes", paste0(gsub("-", "_", gsub("_labeled_count-matrix.csv", "", lf)),
                               "_SN2"))

c1 %>% 
  dplyr::left_join(c2, by = "genes") %>% 
  dplyr::filter(!grepl("FP", genes)) %>% 
  dplyr::arrange(genes) -> dapi_df

dapi_df %>% 
  tidyr::gather(key = "sample", value = "counts", -genes) %>% 
  dplyr::mutate(
    group = dplyr::case_when(
      stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH",
      stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
      stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlTDP43",
      stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
      stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
      stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43"),
    control = dplyr::case_when(
      stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
      stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
      stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
      stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
      stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
      stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43"),
    batch = dplyr::case_when(
      stringr::str_detect(sample, pattern = "SN1") ~ "SN1",
      stringr::str_detect(sample, pattern = "SN2") ~ "SN2"),
    control = as.factor(control),
    control = forcats::fct_relevel(control, c("CtrlFUSH+CtrlTDP43", "FUSH",
                                              "TDP43", "CtrlSOD1", "SOD1"))) %>% 
  dplyr::arrange(genes) -> dapi_df_long

fgl3 <- list()
for (i in 1:length(unique(dapi_df_long$genes))) {
  dapi_df_long %>%
    dplyr::filter(genes == sort(unique(dapi_df_long$genes))[i]) %>%
    ggplot(aes(x = control, y = log(counts), fill = control)) +
    geom_boxplot() +
    geom_point(size = 2, color = "black") +
    ggtitle(paste0(unique(dapi_df_long$genes)[i],
                   " non-normalized pseudobulk counts")) +
    theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1),
          legend.position = "none") -> fgl3[[i]]
}

## save pseudobulk_gene_counts
pdf(file = paste0(fp0, "SN_DAPI_transversal_pseudobulk_gene_count_boxplots.pdf"))
print(fgl3)
dev.off()

################################################################################
########################### NFH transversal SN #################################
################################################################################

############################### SN1 "NFH transversal" #########################
fp1 <- "~/Documents/tmp/2023/202305_piol/count_matrices_correct/SN_Resolve_MC_Images/NFH/non_expanded/"
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
names(c1) <- c("genes", paste0(gsub("-", "_", gsub(
  "_labeled_non-expanded_count-matrix.csv", "", lf)), "_SN1"))

############################### SN2 "NFH transversal" #########################
fp1 <- "~/Documents/tmp/2023/202305_piol/count_matrices_correct/SN_Resolve_Images_2/NFH/"
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
  as.data.frame() -> c2
c2[is.na(c2)] <- 0 ## coerce NA values to zero
names(c2) <- c("genes", paste0(gsub("-", "_", gsub(
  "_labeled_non-expanded_count-matrix.csv", "", lf)), "_SN2"))

c1 %>% 
  dplyr::left_join(c2, by = "genes") %>% 
  dplyr::filter(!grepl("FP", genes)) %>% 
  dplyr::arrange(genes) -> nfh_df

nfh_df %>% 
  tidyr::gather(key = "sample", value = "counts", -genes) %>% 
  dplyr::mutate(
    group = dplyr::case_when(
      stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH",
      stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
      stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlTDP43",
      stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
      stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
      stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43"),
    control = dplyr::case_when(
      stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
      stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
      stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
      stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
      stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
      stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43"),
    batch = dplyr::case_when(
      stringr::str_detect(sample, pattern = "SN1") ~ "SN1",
      stringr::str_detect(sample, pattern = "SN2") ~ "SN2"),
    control = as.factor(control),
    control = forcats::fct_relevel(control, c("CtrlFUSH+CtrlTDP43", "FUSH",
                                              "TDP43", "CtrlSOD1", "SOD1"))) %>% 
  dplyr::arrange(genes) -> nfh_df_long

fgl4 <- list()
for (i in 1:length(unique(nfh_df_long$genes))) {
  nfh_df_long %>%
    dplyr::filter(genes == sort(unique(nfh_df_long$genes))[i]) %>%
    ggplot(aes(x = control, y = log(counts), fill = control)) +
    geom_boxplot() +
    geom_point(size = 2, color = "black") +
    ggtitle(paste0(unique(nfh_df_long$genes)[i],
                   " non-normalized pseudobulk counts")) +
    theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1),
          legend.position = "none") -> fgl4[[i]]
}

pdf(file = paste0(fp0, "SN_NFH_transversal_pseudobulk_gene_count_boxplots.pdf"))
print(fgl4)
dev.off()

################################################################################
########################### CHAT transversal SN ################################
################################################################################

############################### SN1 "CHAT transversal" #########################
fp1 <- "~/Documents/tmp/2023/202305_piol/count_matrices_correct/SN_Resolve_MC_Images/ChAT/non_expanded/"
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
names(c1) <- c("genes", paste0(gsub("-", "_", gsub(
  "_labeled_non-expanded_count-matrix.csv", "", lf)), "_SN1"))

############################### SN2 "CHAT transversal" #########################
fp1 <- "~/Documents/tmp/2023/202305_piol/count_matrices_correct/SN_Resolve_Images_2/ChAT/"
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
  as.data.frame() -> c2
c2[is.na(c2)] <- 0 ## coerce NA values to zero
names(c2) <- c("genes", paste0(gsub("-", "_", gsub(
  "_labeled_non-expanded_count-matrix.csv", "", lf)), "_SN2"))

c1 %>% 
  dplyr::left_join(c2, by = "genes") %>% 
  dplyr::filter(!grepl("FP", genes)) %>% 
  dplyr::arrange(genes) -> chat_df

chat_df %>% 
  tidyr::gather(key = "sample", value = "counts", -genes) %>% 
  dplyr::mutate(
    group = dplyr::case_when(
      stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH",
      stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
      stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlTDP43",
      stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
      stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
      stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43"),
    control = dplyr::case_when(
      stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
      stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
      stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
      stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
      stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
      stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43"),
    batch = dplyr::case_when(
      stringr::str_detect(sample, pattern = "SN1") ~ "SN1",
      stringr::str_detect(sample, pattern = "SN2") ~ "SN2"),
    control = as.factor(control),
    control = forcats::fct_relevel(control, c("CtrlFUSH+CtrlTDP43", "FUSH",
                                              "TDP43", "CtrlSOD1", "SOD1"))) %>% 
  dplyr::arrange(genes) -> chat_df_long

fgl5 <- list()
for (i in 1:length(unique(chat_df_long$genes))) {
  chat_df_long %>%
    dplyr::filter(genes == sort(unique(chat_df_long$genes))[i]) %>%
    ggplot(aes(x = control, y = log(counts), fill = control)) +
    geom_boxplot() +
    geom_point(size = 2, color = "black") +
    ggtitle(paste0(unique(chat_df_long$genes)[i],
                   " non-normalized pseudobulk counts")) +
    theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1),
          legend.position = "none") -> fgl5[[i]]
}

pdf(file = paste0(fp0, "SN_Chat_transversal_pseudobulk_gene_count_boxplots.pdf"))
print(fgl5)
dev.off()

################################################################################
########################### DAPI transversal SC ################################
################################################################################

############################### SN1 "DAPI transversal" #########################
fp1 <- "~/Documents/tmp/2023/202305_piol/count_matrices_correct/SC_Resolve_Images/"
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

c1 %>%
  dplyr::filter(!grepl("FP", genes)) %>% 
  dplyr::arrange(genes) -> dapi_sc_df

dapi_sc_df %>% 
  tidyr::gather(key = "sample", value = "counts", -genes) %>% 
  dplyr::mutate(
    group = dplyr::case_when(
      stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH",
      stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
      stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlTDP43",
      stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
      stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
      stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43"),
    control = dplyr::case_when(
      stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
      stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
      stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
      stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
      stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
      stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43"),
    control = as.factor(control),
    control = forcats::fct_relevel(control, c("CtrlFUSH+CtrlTDP43", "FUSH",
                                              "TDP43", "CtrlSOD1", "SOD1"))) %>% 
  dplyr::arrange(genes) -> dapi_sc_long

fgl6 <- list()
for (i in 1:length(unique(dapi_sc_long$genes))) {
  dapi_sc_long %>%
    dplyr::filter(genes == sort(unique(dapi_sc_long$genes))[i]) %>%
    ggplot(aes(x = control, y = log(counts), fill = control)) +
    geom_boxplot() +
    geom_point(size = 2, color = "black") +
    ggtitle(paste0(unique(dapi_sc_long$genes)[i],
                   " non-normalized pseudobulk counts")) +
    theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1),
          legend.position = "none") -> fgl6[[i]]
}

## save pseudobulk_gene_counts
pdf(file = paste0(fp0, "SC_DAPI_transversal_pseudobulk_gene_count_boxplots.pdf"))
print(fgl6)
dev.off()

################################################################################
########################### DAPI longitudinal SN ###############################
################################################################################

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
names(c1) <- c("genes", paste0(gsub("-", "_", gsub("_labeled_count-matrix.csv", "", lf)),
                               "_SN1"))

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

c1 %>%
  dplyr::left_join(c2, by = "genes") %>%
  dplyr::filter(!grepl("FP", genes)) %>% 
  dplyr::arrange(genes) -> dapi_LONG_df

dapi_LONG_df %>% 
  tidyr::gather(key = "sample", value = "counts", -genes) %>% 
  dplyr::mutate(
    group = dplyr::case_when(
      stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH",
      stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
      stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlTDP43",
      stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
      stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
      stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43"),
    control = dplyr::case_when(
      stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
      stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
      stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
      stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
      stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
      stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43"),
    batch = dplyr::case_when(
      stringr::str_detect(sample, pattern = "SN1") ~ "SN1",
      stringr::str_detect(sample, pattern = "SN2") ~ "SN2"),
    control = as.factor(control),
    control = forcats::fct_relevel(control, c("CtrlFUSH+CtrlTDP43", "FUSH",
                                              "TDP43", "CtrlSOD1", "SOD1"))) %>% 
  dplyr::arrange(genes) -> dapi_df_LONG_long

fgl7 <- list()
for (i in 1:length(unique(dapi_df_LONG_long$genes))) {
  dapi_df_long %>%
    dplyr::filter(genes == sort(unique(dapi_df_LONG_long$genes))[i]) %>%
    ggplot(aes(x = control, y = log(counts), fill = control)) +
    geom_boxplot() +
    geom_point(size = 2, color = "black") +
    ggtitle(paste0(unique(dapi_df_LONG_long$genes)[i],
                   " non-normalized pseudobulk counts")) +
    theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1),
          legend.position = "none") -> fgl7[[i]]
}

## save pseudobulk_gene_counts
pdf(file = paste0(fp0, "SN_DAPI_longitudinal_pseudobulk_gene_count_boxplots.pdf"))
print(fgl7)
dev.off()

############################## save pseudobulk counts ##########################
count_list <- list(dapi_df, dapi_LONG_df, nfh_df, chat_df,
                   dapi_sc_df, NFH_1500_total_counts, new_seg_df)
names(count_list) <- c("DAPI SN transversal", "DAPI SN longitudinal",
                       "NFH SN transversal", "Chat SN transversal",
                       "DAPI SC transversal",
                       "SN NFH 1500 Chatpos vs Chatneg",
                       "NFH vs segmentation")
writexl::write_xlsx(x = count_list, path = paste0(fp0, "all_pseudobulk_counts.xlsx"))
