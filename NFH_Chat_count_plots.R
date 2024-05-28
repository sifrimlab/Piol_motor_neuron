## NFH 1500 gene count plots

## load libraries
library("dplyr")
library("ggplot2")
library("readr")
library("writexl")
library("tibble")
library("stringr")
library("forcats")

## there the ggplot PDF will be saved
out <- "~/Documents/tmp/dacruz/202311_piol/"

############################# total counts #####################################

######################### read "SN1 Chat_trans"  
# file path input files
fp <- paste0("~/Documents/tmp/dacruz/202305_piol/count_matrices_correct/with_area/",
             "SN_Resolve_MC_Images/non_expanded/ChAT/")
lf <- list.files(path = fp, pattern = "*.csv")
dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    tibble::column_to_rownames(var = "cell_label") %>%
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
  "_labeled_non-expanded_count-matrix_with-area.csv", "", lf)), "_SN1"))

######################### read "SN2 Chat_trans"
## file path input files
fp <- paste0("~/Documents/tmp/dacruz/202305_piol/count_matrices_correct/with_area/",
             "SN_Resolve_Images_2/non_expanded/ChAT/")
lf <- list.files(path = fp, pattern = "*.csv")
dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    tibble::column_to_rownames(var = "cell_label") %>%
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
  "_labeled_non-expanded_count-matrix_with-area.csv", "", lf)), "_SN2"))

cp1 %>%
  dplyr::left_join(cp2, by = "genes") %>%
  tidyr::gather(key = "sample", value = "counts", -genes) -> total_counts

total_counts %>%
  dplyr::group_by(genes) %>%
  dplyr::summarise(mean_counts = mean(counts, na.rm = TRUE)) -> mean_counts

total_counts %>%
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
  dplyr::left_join(mean_counts, by = "genes") %>%
  dplyr::arrange(genes) -> all_Chat

fgl1 <- list()
for (i in 1:length(unique(all_Chat$genes))) {
  all_Chat %>%
    dplyr::filter(genes == unique(all_Chat$genes)[i]) %>%
    ggplot(aes(x = control, y = counts, color = control)) +
    geom_point() +
    geom_hline(yintercept = all_Chat %>%
                 dplyr::filter(genes == unique(all_Chat$genes)[i]) %>%
                 dplyr::slice(1) %>%
                 dplyr::pull(mean_counts),
               color = "red", linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1),
          legend.position = "none") +
    ggtitle(paste0(unique(all_Chat$genes)[i],
                   " non-normalized pseudobulk counts in SN Chat")) -> fgl1[[i]]
}

# save pseudobulk_gene_counts
pdf(file = paste0(out, "SN_Chat_total_gene_counts.pdf"))
print(fgl1)
dev.off()

######################### read "SN1 NFH_trans"  
# file path input files
# fp <- paste0("~/Documents/tmp/dacruz/202305_piol/count_matrices_correct/with_area/",
#              "SN_Resolve_MC_Images/non_expanded/NFH/")
# lf <- list.files(path = fp, pattern = "*.csv")
# dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
# dfs2 <- lapply(dfs1, function(i){i %>%
#     tibble::column_to_rownames(var = "cell_label") %>%
#     # dplyr::filter(area > 1500) %>%
#     dplyr::select(-c(area, contains("FP")))})
# dfs3 <- lapply(dfs2, function(i){t(i)})
# rl1 <- list()
# for (i in 1:length(dfs3)) {
#   rl1[[i]] <- as.data.frame(rowSums(dfs3[[i]]))
# }
# tibble::rownames_to_column(rl1[[1]], var = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[2]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[3]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[4]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[5]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[6]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[7]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[8]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[9]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[10]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[11]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[12]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[13]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[14]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[15]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[16]], var = "genes"), by = "genes") %>%
#   as.data.frame() -> cp1
# cp1[is.na(cp1)] <- 0 ## coerce NA values to zero
# names(cp1) <- c("genes", paste0(gsub("-", "_", gsub(
#   "_labeled_non-expanded_count-matrix_with-area.csv", "", lf)), "_SN1"))
# 
# ######################### read "SN2 NFH_trans"
# ## file path input files
# fp <- paste0("~/Documents/tmp/dacruz/202305_piol/count_matrices_correct/with_area/",
#              "SN_Resolve_Images_2/non_expanded/NFH/")
# lf <- list.files(path = fp, pattern = "*.csv")
# dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
# dfs2 <- lapply(dfs1, function(i){i %>%
#     tibble::column_to_rownames(var = "cell_label") %>%
#     # dplyr::filter(area > 1500) %>%
#     dplyr::select(-c(area, contains("FP")))})
# dfs3 <- lapply(dfs2, function(i){t(i)})
# rl1 <- list()
# for (i in 1:length(dfs3)) {
#   rl1[[i]] <- as.data.frame(rowSums(dfs3[[i]]))
# }
# tibble::rownames_to_column(rl1[[1]], var = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[2]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[3]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[4]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[5]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[6]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[7]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[8]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[9]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[10]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[11]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[12]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[13]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[14]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[15]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl1[[16]], var = "genes"), by = "genes") %>%
#   as.data.frame() -> cp2
# cp2[is.na(cp2)] <- 0 ## coerce NA values to zero
# names(cp2) <- c("genes", paste0(gsub("-", "_", gsub(
#   "_labeled_non-expanded_count-matrix_with-area.csv", "", lf)), "_SN2"))
# 
# cp1 %>%
#   dplyr::left_join(cp2, by = "genes") %>%
#   tidyr::gather(key = "sample", value = "counts", -genes) -> total_counts
# 
# total_counts %>%
#     dplyr::group_by(genes) %>%
#     dplyr::summarise(mean_counts = mean(counts, na.rm = TRUE)) -> mean_counts
# 
# total_counts %>%
#   dplyr::mutate(
#     group = dplyr::case_when(
#       stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH",
#       stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
#       stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlTDP43",
#       stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
#       stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
#       stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43"),
#     control = dplyr::case_when(
#       stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
#       stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
#       stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
#       stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
#       stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
#       stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43"),
#     batch = dplyr::case_when(
#       stringr::str_detect(sample, pattern = "SN1") ~ "SN1",
#       stringr::str_detect(sample, pattern = "SN2") ~ "SN2"),
#     control = as.factor(control),
#     control = forcats::fct_relevel(control, c("CtrlFUSH+CtrlTDP43", "FUSH",
#                                               "TDP43", "CtrlSOD1", "SOD1"))) %>%
#   dplyr::left_join(mean_counts, by = "genes") %>%
#   dplyr::arrange(genes) -> all_NFH
# 
# fgl1 <- list()
# for (i in 1:length(unique(all_NFH$genes))) {
#   all_NFH %>%
#     dplyr::filter(genes == unique(all_NFH$genes)[i]) %>%
#     ggplot(aes(x = control, y = counts, color = control)) +
#     geom_point() +
#     geom_hline(yintercept = all_NFH %>%
#                  dplyr::filter(genes == unique(all_NFH$genes)[i]) %>%
#                  dplyr::slice(1) %>%
#                  dplyr::pull(mean_counts),
#                color = "red", linetype = "dashed") +
#     theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1),
#           legend.position = "none") +
#     ggtitle(paste0(unique(all_NFH$genes)[i],
#                    " non-normalized pseudobulk counts in SN NFH")) -> fgl1[[i]]
# }
# 
# # save pseudobulk_gene_counts
# pdf(file = paste0(out, "SN_NFH_total_gene_counts.pdf"))
# print(fgl1)
# dev.off()

############################# segmented counts ##################################

######################### read "SN1 NFH_trans"  
# file path input files
# fp <- paste0("~/Documents/tmp/dacruz/202305_piol/count_matrices_correct/with_area/",
#              "SN_Resolve_MC_Images/non_expanded/NFH/")
# lf <- list.files(path = fp, pattern = "*.csv")
# dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
# dfs2 <- lapply(dfs1, function(i){i %>%
#     tibble::column_to_rownames(var = "cell_label") %>%
#     dplyr::select(-c(area, contains("FP"))) %>% 
#     tibble::rownames_to_column(var = "cell_id") %>% 
#     dplyr::mutate(cell_id = paste0(
#       gsub("-", "_",
#           gsub("_labeled_non-expanded_count-matrix_with-area.csv", "", lf)),
#               "_", cell_id, "_SN1")) %>% 
#     tibble::column_to_rownames(var = "cell_id")
#   })
# dfs3 <- lapply(dfs2, function(i){t(i) %>%
#     as.data.frame() %>%
#     tibble::rownames_to_column(var = "genes")})
# dfs3[[1]] %>%
#   dplyr::full_join(dfs3[[2]], by = "genes") %>%
#   dplyr::full_join(dfs3[[3]], by = "genes") %>%
#   dplyr::full_join(dfs3[[4]], by = "genes") %>%
#   dplyr::full_join(dfs3[[5]], by = "genes") %>%
#   dplyr::full_join(dfs3[[6]], by = "genes") %>%
#   dplyr::full_join(dfs3[[7]], by = "genes") %>%
#   dplyr::full_join(dfs3[[8]], by = "genes") %>%
#   dplyr::full_join(dfs3[[9]], by = "genes") %>%
#   dplyr::full_join(dfs3[[10]], by = "genes") %>%
#   dplyr::full_join(dfs3[[11]], by = "genes") %>%
#   dplyr::full_join(dfs3[[12]], by = "genes") %>%
#   dplyr::full_join(dfs3[[13]], by = "genes") %>%
#   dplyr::full_join(dfs3[[14]], by = "genes") %>%
#   dplyr::full_join(dfs3[[15]], by = "genes") %>%
#   dplyr::full_join(dfs3[[16]], by = "genes") %>%
#   as.data.frame() -> cv1
# 
# ######################### read "SN2 NFH_trans"
# ## file path input files
# fp <- paste0("~/Documents/tmp/dacruz/202305_piol/count_matrices_correct/with_area/",
#              "SN_Resolve_Images_2/non_expanded/NFH/")
# lf <- list.files(path = fp, pattern = "*.csv")
# dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
# dfs2 <- lapply(dfs1, function(i){i %>%
#     tibble::column_to_rownames(var = "cell_label") %>%
#     dplyr::select(-c(area, contains("FP"))) %>% 
#     tibble::rownames_to_column(var = "cell_id") %>% 
#     dplyr::mutate(cell_id = paste0(
#       gsub("-", "_",
#            gsub("_labeled_non-expanded_count-matrix_with-area.csv", "", lf)),
#       "_", cell_id, "_SN2")) %>% 
#     tibble::column_to_rownames(var = "cell_id")
# })
# dfs3 <- lapply(dfs2, function(i){t(i) %>%
#     as.data.frame() %>%
#     tibble::rownames_to_column(var = "genes")})
# dfs3[[1]] %>%
#   dplyr::full_join(dfs3[[2]], by = "genes") %>%
#   dplyr::full_join(dfs3[[3]], by = "genes") %>%
#   dplyr::full_join(dfs3[[4]], by = "genes") %>%
#   dplyr::full_join(dfs3[[5]], by = "genes") %>%
#   dplyr::full_join(dfs3[[6]], by = "genes") %>%
#   dplyr::full_join(dfs3[[7]], by = "genes") %>%
#   dplyr::full_join(dfs3[[8]], by = "genes") %>%
#   dplyr::full_join(dfs3[[9]], by = "genes") %>%
#   dplyr::full_join(dfs3[[10]], by = "genes") %>%
#   dplyr::full_join(dfs3[[11]], by = "genes") %>%
#   dplyr::full_join(dfs3[[12]], by = "genes") %>%
#   dplyr::full_join(dfs3[[13]], by = "genes") %>%
#   dplyr::full_join(dfs3[[14]], by = "genes") %>%
#   dplyr::full_join(dfs3[[15]], by = "genes") %>%
#   dplyr::full_join(dfs3[[16]], by = "genes") %>%
#   as.data.frame() -> cv2
# 
# cv1 %>%
#   dplyr::left_join(cv2, by = "genes") %>%
#   tidyr::gather(key = "sample", value = "counts", -genes) -> total_counts
# 
# total_counts %>%
#     dplyr::group_by(genes) %>%
#     dplyr::summarise(mean_counts = mean(counts, na.rm = TRUE)) -> mean_counts
# 
# total_counts %>%
#   dplyr::mutate(
#     group = dplyr::case_when(
#       stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH",
#       stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
#       stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlTDP43",
#       stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
#       stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
#       stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43"),
#     control = dplyr::case_when(
#       stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
#       stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
#       stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
#       stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
#       stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
#       stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43"),
#     batch = dplyr::case_when(
#       stringr::str_detect(sample, pattern = "SN1") ~ "SN1",
#       stringr::str_detect(sample, pattern = "SN2") ~ "SN2"),
#     control = as.factor(control),
#     control = forcats::fct_relevel(control, c("CtrlFUSH+CtrlTDP43", "FUSH",
#                                               "TDP43", "CtrlSOD1", "SOD1"))) %>%
#   dplyr::left_join(mean_counts, by = "genes") %>%
#   dplyr::arrange(genes) -> all_NFH
# 
# fgl2 <- list()
# for (i in 1:length(unique(all_NFH$genes))) {
#   all_NFH %>%
#     dplyr::filter(genes == unique(all_NFH$genes)[i]) %>%
#     ggplot(aes(x = control, y = counts, color = control)) +
#     geom_jitter() +
#     geom_hline(yintercept = all_NFH %>%
#                  dplyr::filter(genes == unique(all_NFH$genes)[i]) %>%
#                  dplyr::slice(1) %>%
#                  dplyr::pull(mean_counts),
#                color = "red", linetype = "dashed") +
#     theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1),
#           legend.position = "none") +
#     ggtitle(paste0(unique(all_NFH$genes)[i],
#                    " non-normalized segmented counts in SN NFH")) -> fgl2[[i]]
# }
# 
# # save pseudobulk_gene_counts
# pdf(file = paste0(out, "SN_NFH_segmented_gene_counts.pdf"))
# print(fgl2)
# dev.off()