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

## segments with an anomylous count value 
## CtrlFUSH_1_A2_2_NFH lf3
## CtrlFUSH_1_A2_2_NFH_row3221_col3215_7786_SN1
## CtrlFUSH_1_A2_2_NFH_row3259_col3212_1125_SN1
## CtrlFUSH_1_A2_2_ChAT_row3221_col3215_5360_SN1
## CtrlFUSH_1_A2_2_ChAT_row3232_col3213_27559_SN1

######################### NFH Chatpos ( area > 1500 ) SN1
fp1 <- paste0("~/Documents/tmp/dacruz/202305_piol/count_matrices_correct/with_area/",
             "SN_Resolve_MC_Images/non_expanded/NFH/")
lf1 <- list.files(path = fp1, pattern = "*.csv")
dfs1 <- lapply(lf1, function(i){readr::read_csv(file = paste0(fp1, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    dplyr::filter(!grepl("row3221_col3215", cell_label),
                  !grepl("row3259_col3212", cell_label),
                  !grepl("row3221_col3215", cell_label),
                  !grepl("row3232_col3213", cell_label)) %>%
    dplyr::filter(area > 1500) %>%
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
  as.data.frame() -> chat1
chat1[is.na(chat1)] <- 0 ## coerce NA values to zero
names(chat1) <- c("genes", paste0(gsub("-", "_", gsub(
  "_labeled_non-expanded_count-matrix_with-area.csv", "", lf1)), "_SN1_Chatpos"))

######################### NFH Chatpos ( area > 1500 ) SN2
fp2 <- paste0("~/Documents/tmp/dacruz/202305_piol/count_matrices_correct/with_area/",
             "SN_Resolve_Images_2/non_expanded/NFH/")
lf2 <- list.files(path = fp2, pattern = "*.csv")
dfs1 <- lapply(lf2, function(i){readr::read_csv(file = paste0(fp2, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    dplyr::filter(!grepl("row3221_col3215", cell_label),
                  !grepl("row3259_col3212", cell_label),
                  !grepl("row3221_col3215", cell_label),
                  !grepl("row3232_col3213", cell_label)) %>%
    dplyr::filter(area > 1500) %>%
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
  as.data.frame() -> chat2
chat2[is.na(chat2)] <- 0 ## coerce NA values to zero
names(chat2) <- c("genes", paste0(gsub("-", "_", gsub(
  "_labeled_non-expanded_count-matrix_with-area.csv", "", lf2)), "_SN2_Chatpos"))

chat1 %>%
  dplyr::left_join(chat2, by = "genes") %>%
  tibble::column_to_rownames(var = "genes") %>%
  as.data.frame() -> total_counts

data.frame(colSums(total_counts)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  dplyr::rename(all_counts = names(.)[2]) %>%
  dplyr::mutate(
    motor_neuron = dplyr::case_when(
      grepl("Chatpos", sample) ~ "Chat_pos",
      grepl("Chatneg", sample) ~ "Chat_neg"),
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
                                              "TDP43", "CtrlSOD1", "SOD1"))
    ) -> all_chatpos

######################### NFH Chatneg ( area < 1500 ) SN1
fp3 <- paste0("~/Documents/tmp/dacruz/202305_piol/count_matrices_correct/with_area/",
             "SN_Resolve_MC_Images/non_expanded/NFH/")
lf3 <- list.files(path = fp3, pattern = "*.csv")
dfs1 <- lapply(lf3, function(i){readr::read_csv(file = paste0(fp3, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    dplyr::filter(!grepl("row3221_col3215", cell_label),
                  !grepl("row3259_col3212", cell_label),
                  !grepl("row3221_col3215", cell_label),
                  !grepl("row3232_col3213", cell_label)) %>%
    dplyr::filter(area < 1500) %>%
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
  as.data.frame() -> nfh1
nfh1[is.na(nfh1)] <- 0 ## coerce NA values to zero
names(nfh1) <- c("genes", paste0(gsub("-", "_", gsub(
  "_labeled_non-expanded_count-matrix_with-area.csv", "", lf3)), "_SN1_Chatneg"))

######################### NFH Chatneg ( area < 1500 ) SN2
fp4 <- paste0("~/Documents/tmp/dacruz/202305_piol/count_matrices_correct/with_area/",
             "SN_Resolve_Images_2/non_expanded/NFH/")
lf4 <- list.files(path = fp4, pattern = "*.csv")
dfs1 <- lapply(lf4, function(i){readr::read_csv(file = paste0(fp4, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    dplyr::filter(!grepl("row3221_col3215", cell_label),
                  !grepl("row3259_col3212", cell_label),
                  !grepl("row3221_col3215", cell_label),
                  !grepl("row3232_col3213", cell_label)) %>%
    dplyr::filter(area < 1500) %>%
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
  as.data.frame() -> nfh2
nfh2[is.na(nfh2)] <- 0 ## coerce NA values to zero
names(nfh2) <- c("genes", paste0(gsub("-", "_", gsub(
  "_labeled_non-expanded_count-matrix_with-area.csv", "", lf4)), "_SN2_Chatneg"))

nfh1 %>%
  dplyr::left_join(nfh2, by = "genes") %>%
  tibble::column_to_rownames(var = "genes") %>%
  as.data.frame() -> total_counts

data.frame(colSums(total_counts)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  dplyr::rename(all_counts = names(.)[2]) %>%
  dplyr::mutate(
    motor_neuron = dplyr::case_when(
      grepl("Chatpos", sample) ~ "Chat_pos",
      grepl("Chatneg", sample) ~ "Chat_neg"),
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
                                              "TDP43", "CtrlSOD1", "SOD1"))
  ) -> all_chatneg

dplyr::bind_rows(all_chatpos, all_chatneg) -> all_non_seg
all_non_seg %>%
  ggplot(aes(x = control, y = all_counts, color = control)) +
  geom_point() +
  geom_hline(yintercept = mean(all_chatpos$all_counts, na.rm = TRUE),
             color = "red", linetype = "dashed") +
  geom_hline(yintercept = mean(all_chatneg$all_counts, na.rm = TRUE),
             color = "blue", linetype = "dashed") +
  theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1),
        legend.position = "none") +
  ylab("Total gene counts per sample") +
  ggtitle("Total combined non-normalized pseudobulk NFH counts in sciatic nerve") +
  facet_grid(~motor_neuron) -> p2a

############################# segmented counts #################################

######################### NFH Chatpos ( area > 1500 ) SN1
fp1 <- paste0("~/Documents/tmp/dacruz/202305_piol/count_matrices_correct/with_area/",
             "SN_Resolve_MC_Images/non_expanded/NFH/")
lf1 <- list.files(path = fp1, pattern = "*.csv")
dfs1 <- lapply(lf1, function(i){readr::read_csv(file = paste0(fp1, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    tibble::column_to_rownames(var = "cell_label") %>%
    dplyr::filter(area > 1500) %>%
    dplyr::select(-c(area, contains("FP"))) %>%
    tibble::rownames_to_column(var = "cell_id") %>%
    dplyr::mutate(cell_id = paste0(
      gsub("-", "_",
           gsub("_labeled_non-expanded_count-matrix_with-area.csv", "", lf1)),
      "_", cell_id, "_SN1_Chatpos")) %>%
    tibble::column_to_rownames(var = "cell_id")
})
dfs3 <- lapply(dfs2, function(i){t(i) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "genes")})
dfs3[[1]] %>%
  dplyr::full_join(dfs3[[2]], by = "genes") %>%
  dplyr::full_join(dfs3[[3]], by = "genes") %>%
  dplyr::full_join(dfs3[[4]], by = "genes") %>%
  dplyr::full_join(dfs3[[5]], by = "genes") %>%
  dplyr::full_join(dfs3[[6]], by = "genes") %>%
  dplyr::full_join(dfs3[[7]], by = "genes") %>%
  dplyr::full_join(dfs3[[8]], by = "genes") %>%
  dplyr::full_join(dfs3[[9]], by = "genes") %>%
  dplyr::full_join(dfs3[[10]], by = "genes") %>%
  dplyr::full_join(dfs3[[11]], by = "genes") %>%
  dplyr::full_join(dfs3[[12]], by = "genes") %>%
  dplyr::full_join(dfs3[[13]], by = "genes") %>%
  dplyr::full_join(dfs3[[14]], by = "genes") %>%
  dplyr::full_join(dfs3[[15]], by = "genes") %>%
  dplyr::full_join(dfs3[[16]], by = "genes") %>%
  as.data.frame() -> chat1s

######################### NFH Chatpos ( area > 1500 ) SN2
fp2 <- paste0("~/Documents/tmp/dacruz/202305_piol/count_matrices_correct/with_area/",
             "SN_Resolve_Images_2/non_expanded/NFH/")
lf2 <- list.files(path = fp2, pattern = "*.csv")
dfs1 <- lapply(lf2, function(i){readr::read_csv(file = paste0(fp2, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    tibble::column_to_rownames(var = "cell_label") %>%
    dplyr::filter(area > 1500) %>%
    dplyr::select(-c(area, contains("FP"))) %>%
    tibble::rownames_to_column(var = "cell_id") %>%
    dplyr::mutate(cell_id = paste0(
      gsub("-", "_",
           gsub("_labeled_non-expanded_count-matrix_with-area.csv", "", lf2)),
      "_", cell_id, "_SN2_Chatpos")) %>%
    tibble::column_to_rownames(var = "cell_id")
})
dfs3 <- lapply(dfs2, function(i){t(i) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "genes")})
dfs3[[1]] %>%
  dplyr::full_join(dfs3[[2]], by = "genes") %>%
  dplyr::full_join(dfs3[[3]], by = "genes") %>%
  dplyr::full_join(dfs3[[4]], by = "genes") %>%
  dplyr::full_join(dfs3[[5]], by = "genes") %>%
  dplyr::full_join(dfs3[[6]], by = "genes") %>%
  dplyr::full_join(dfs3[[7]], by = "genes") %>%
  dplyr::full_join(dfs3[[8]], by = "genes") %>%
  dplyr::full_join(dfs3[[9]], by = "genes") %>%
  dplyr::full_join(dfs3[[10]], by = "genes") %>%
  dplyr::full_join(dfs3[[11]], by = "genes") %>%
  dplyr::full_join(dfs3[[12]], by = "genes") %>%
  dplyr::full_join(dfs3[[13]], by = "genes") %>%
  dplyr::full_join(dfs3[[14]], by = "genes") %>%
  dplyr::full_join(dfs3[[15]], by = "genes") %>%
  dplyr::full_join(dfs3[[16]], by = "genes") %>%
  as.data.frame() -> chat2s

chat1s %>%
  dplyr::left_join(chat2s, by = "genes") %>%
  tibble::column_to_rownames(var = "genes") %>%
  as.data.frame() -> total_counts

total_counts[is.na(total_counts)] <- 0

data.frame(colSums(total_counts)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  dplyr::rename(all_counts = names(.)[2]) %>%
  dplyr::mutate(
    motor_neuron = dplyr::case_when(
      grepl("Chatpos", sample) ~ "Chat_pos",
      grepl("Chatneg", sample) ~ "Chat_neg"),
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
      dplyr::filter(all_counts < 300) -> all_chatpos_seg

######################### read NFH Chatneg ( area < 1500 ) SN1
fp3 <- paste0("~/Documents/tmp/dacruz/202305_piol/count_matrices_correct/with_area/",
              "SN_Resolve_MC_Images/non_expanded/NFH/")
lf3 <- list.files(path = fp3, pattern = "*.csv")
dfs1 <- lapply(lf3, function(i){readr::read_csv(file = paste0(fp3, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    tibble::column_to_rownames(var = "cell_label") %>%
    dplyr::filter(area < 1500) %>%
    dplyr::select(-c(area, contains("FP"))) %>%
    tibble::rownames_to_column(var = "cell_id") %>%
    dplyr::mutate(cell_id = paste0(
      gsub("-", "_",
           gsub("_labeled_non-expanded_count-matrix_with-area.csv", "", lf3)),
      "_", cell_id, "_SN1_Chatneg")) %>%
    tibble::column_to_rownames(var = "cell_id")
})
dfs3 <- lapply(dfs2, function(i){t(i) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "genes")})
dfs3[[1]] %>%
  dplyr::full_join(dfs3[[2]], by = "genes") %>%
  dplyr::full_join(dfs3[[3]], by = "genes") %>%
  dplyr::full_join(dfs3[[4]], by = "genes") %>%
  dplyr::full_join(dfs3[[5]], by = "genes") %>%
  dplyr::full_join(dfs3[[6]], by = "genes") %>%
  dplyr::full_join(dfs3[[7]], by = "genes") %>%
  dplyr::full_join(dfs3[[8]], by = "genes") %>%
  dplyr::full_join(dfs3[[9]], by = "genes") %>%
  dplyr::full_join(dfs3[[10]], by = "genes") %>%
  dplyr::full_join(dfs3[[11]], by = "genes") %>%
  dplyr::full_join(dfs3[[12]], by = "genes") %>%
  dplyr::full_join(dfs3[[13]], by = "genes") %>%
  dplyr::full_join(dfs3[[14]], by = "genes") %>%
  dplyr::full_join(dfs3[[15]], by = "genes") %>%
  dplyr::full_join(dfs3[[16]], by = "genes") %>%
  as.data.frame() -> nfh1s
# nfh1s[is.na(nfh1s)] <- 0 ## coerce NA values to zero

######################### read NFH Chatneg ( area < 1500 ) SN2
fp4 <- paste0("~/Documents/tmp/dacruz/202305_piol/count_matrices_correct/with_area/",
              "SN_Resolve_Images_2/non_expanded/NFH/")
lf4 <- list.files(path = fp4, pattern = "*.csv")
dfs1 <- lapply(lf4, function(i){readr::read_csv(file = paste0(fp4, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    tibble::column_to_rownames(var = "cell_label") %>%
    dplyr::filter(area < 1500) %>%
    dplyr::select(-c(area, contains("FP"))) %>%
    tibble::rownames_to_column(var = "cell_id") %>%
    dplyr::mutate(cell_id = paste0(
      gsub("-", "_",
           gsub("_labeled_non-expanded_count-matrix_with-area.csv", "", lf4)),
      "_", cell_id, "_SN2_Chatneg")) %>%
    tibble::column_to_rownames(var = "cell_id")
})
dfs3 <- lapply(dfs2, function(i){t(i) %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "genes")})
dfs3[[1]] %>%
  dplyr::full_join(dfs3[[2]], by = "genes") %>%
  dplyr::full_join(dfs3[[3]], by = "genes") %>%
  dplyr::full_join(dfs3[[4]], by = "genes") %>%
  dplyr::full_join(dfs3[[5]], by = "genes") %>%
  dplyr::full_join(dfs3[[6]], by = "genes") %>%
  dplyr::full_join(dfs3[[7]], by = "genes") %>%
  dplyr::full_join(dfs3[[8]], by = "genes") %>%
  dplyr::full_join(dfs3[[9]], by = "genes") %>%
  dplyr::full_join(dfs3[[10]], by = "genes") %>%
  dplyr::full_join(dfs3[[11]], by = "genes") %>%
  dplyr::full_join(dfs3[[12]], by = "genes") %>%
  dplyr::full_join(dfs3[[13]], by = "genes") %>%
  dplyr::full_join(dfs3[[14]], by = "genes") %>%
  dplyr::full_join(dfs3[[15]], by = "genes") %>%
  dplyr::full_join(dfs3[[16]], by = "genes") %>%
  as.data.frame() -> nfh2s
# nfh2s[is.na(nfh2s)] <- 0 ## coerce NA values to zero

nfh1s %>%
  dplyr::left_join(nfh2s, by = "genes") %>%
  tibble::column_to_rownames(var = "genes") %>%
  as.data.frame() -> total_counts

total_counts[is.na(total_counts)] <- 0

data.frame(colSums(total_counts)) %>%
  tibble::rownames_to_column(var = "sample") %>%
  dplyr::rename(all_counts = names(.)[2]) %>%
  dplyr::mutate(
    motor_neuron = dplyr::case_when(
      grepl("Chatpos", sample) ~ "Chat_pos",
      grepl("Chatneg", sample) ~ "Chat_neg"),
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
  dplyr::filter(all_counts < 300) -> all_chatneg_seg

dplyr::bind_rows(all_chatpos_seg, all_chatneg_seg) -> all_seg
all_seg %>% 
    ggplot(aes(x = control, y = log10(all_counts), color = control, fill = control)) +
    geom_jitter(size = 0.5, alpha = 0.5, color = "black") +
    geom_violin() +
    geom_hline(yintercept = log10(mean(all_chatpos_seg$all_counts, na.rm = TRUE)),
               color = "red", linetype = "dashed") +
    geom_hline(yintercept = log10(mean(all_chatneg_seg$all_counts, na.rm = TRUE)),
               color = "blue", linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1),
          legend.position = "none") +
    ylab("log10 total gene counts per segment") +
    ggtitle("Total non-normalized pseudobulk counts per NFH segment in sciatic nerve") +
    facet_grid(~motor_neuron) -> p4a

## save pseudobulk_gene_counts
pdf(file = paste0(out, "SN_NFH_1500_segment_gene_counts.pdf"))
p2a; p4a
dev.off()
