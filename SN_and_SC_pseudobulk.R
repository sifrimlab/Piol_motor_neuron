## This script aggregates cell segmented count data from ALS Resolve data for
## DAPI, ChAT and NFH stains into pseudobulk count matrices by full joining row
## sum counts of each raw cell segmentation count matrix. These pseudo matrices
## are collated and analyzes to produce various QC visualizations. The processed
## counts are saved in Excel format

## references
## https://stackoverflow.com/questions/40315227/how-to-solve-prcomp-default-cannot-rescale-a-constant-zero-column-to-unit-var#:~:text=The%20error%20is%20because%20one,find%20the%20zero%20variance%20variables.
## https://bioconductor.org/packages/devel/bioc/vignettes/PCAtools/inst/doc/PCAtools.html

## load libraries
library("dplyr")
library("readr")
library("writexl")
library("tibble")
library("stringr")
library("forcats")
library("ggplot2")
library("ggrepel")
library("ggfortify")

## other parameters
options(scipen = 999)
fp0 <- "~/Documents/tmp/2023/202307_piol/data/RDS/"

##################################### SC ######################################
################################# "SC DAPI_trans" #############################

fp <- paste0("~/Documents/tmp/2023/202305_piol/count_matrices_correct/with_area/",
             "SC_Resolve_Images/")
lf <- list.files(path = fp, pattern = "*.csv")
dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    dplyr::filter(area > 2000,
                  area < 15000) %>%
    dplyr::mutate(total_objects = 1) %>%
    dplyr::rename(total_area = area) %>%
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
  as.data.frame() -> c1
c1[is.na(c1)] <- 0 ## coerce NA values to zero
names(c1) <- c("genes", paste0(
  gsub("-", "_", gsub("_labeled_count-matrix_with-area.csv", "", lf)), "_SC"))

##################################### SN1 ######################################
################################# "SN1 DAPI_trans" #############################
fp <- paste0("~/Documents/tmp/2023/202305_piol/count_matrices_correct/with_area/",
             "SN_Resolve_MC_Images/expanded/DAPI/transversal/")
lf <- list.files(path = fp, pattern = "*.csv")
dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    dplyr::filter(area > 2000,
                  area < 15000) %>%
    dplyr::mutate(total_objects = 1) %>%
    dplyr::rename(total_area = area) %>%
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
  as.data.frame() -> c2a
c2a[is.na(c2a)] <- 0 ## coerce NA values to zero
names(c2a) <- c("genes", paste0(
  gsub("-", "_", gsub("_labeled_count-matrix_with-area.csv", "", lf)), "_SN1"))

############################ "SN1 ChAT_transversal" ############################
fp <- paste0("~/Documents/tmp/2023/202305_piol/count_matrices_correct/with_area/",
             "SN_Resolve_MC_Images/non_expanded/ChAT/")
lf <- list.files(path = fp, pattern = "*.csv")
dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    dplyr::filter(area < 2500) %>%
    dplyr::mutate(total_objects = 1) %>%
    dplyr::rename(total_area = area) %>%
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
  as.data.frame() -> c3a
c3a[is.na(c3a)] <- 0 ## coerce NA values to zero
names(c3a) <- c("genes", paste0(
  gsub("-", "_", gsub("_labeled_non-expanded_count-matrix_with-area.csv", "", lf)), "_SN1"))

############################ "SN1 NFH_transversal" ############################
fp <- paste0("~/Documents/tmp/2023/202305_piol/count_matrices_correct/with_area/",
             "SN_Resolve_MC_Images/non_expanded/NFH/")
lf <- list.files(path = fp, pattern = "*.csv")
dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    dplyr::filter(area < 2500) %>%
    dplyr::mutate(total_objects = 1) %>%
    dplyr::rename(total_area = area) %>%
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
  as.data.frame() -> c4a
c4a[is.na(c4a)] <- 0 ## coerce NA values to zero
names(c4a) <- c("genes", paste0(
  gsub("-", "_", gsub("_labeled_non-expanded_count-matrix_with-area.csv", "", lf)), "_SN1"))

#################################### SN2 #######################################
################################# "SN2 DAPI_trans" #############################
fp <- paste0("~/Documents/tmp/2023/202305_piol/count_matrices_correct/with_area/",
             "SN_Resolve_Images_2/expanded/DAPI_transversal/")
lf <- list.files(path = fp, pattern = "*.csv")
dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    dplyr::filter(area > 2000,
                  area < 15000) %>%
    dplyr::mutate(total_objects = 1) %>%
    dplyr::rename(total_area = area) %>%
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
  as.data.frame() -> c2b
c2b[is.na(c2b)] <- 0 ## coerce NA values to zero
names(c2b) <- c("genes", paste0(
  gsub("-", "_", gsub("_labeled_count-matrix_with-area.csv", "", lf)), "_SN2"))

############################ "SN2 ChAT_transversal" ############################
fp <- paste0("~/Documents/tmp/2023/202305_piol/count_matrices_correct/with_area/",
             "SN_Resolve_Images_2/non_expanded/ChAT/")
lf <- list.files(path = fp, pattern = "*.csv")
dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    dplyr::filter(area < 2500) %>%
    dplyr::mutate(total_objects = 1) %>%
    dplyr::rename(total_area = area) %>%
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
  as.data.frame() -> c3b
c3b[is.na(c3b)] <- 0 ## coerce NA values to zero
names(c3b) <- c("genes", paste0(
  gsub("-", "_", gsub("_labeled_non-expanded_count-matrix_with-area.csv", "", lf)), "_SN2"))

############################# "SN2 NFH_transversal" ############################
fp <- paste0("~/Documents/tmp/2023/202305_piol/count_matrices_correct/with_area/",
             "SN_Resolve_Images_2/non_expanded/NFH/")
lf <- list.files(path = fp, pattern = "*.csv")
dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
dfs2 <- lapply(dfs1, function(i){i %>%
    dplyr::filter(area < 2500) %>%
    dplyr::mutate(total_objects = 1) %>%
    dplyr::rename(total_area = area) %>%
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
  as.data.frame() -> c4b
c4b[is.na(c4b)] <- 0 ## coerce NA values to zero
names(c4b) <- c("genes", paste0(
  gsub("-", "_", gsub("_labeled_non-expanded_count-matrix_with-area.csv", "", lf)), "_SN2"))

################################# COLLATE THE MATRICES #########################
c1 %>%
  dplyr::mutate(total_counts = rowSums(across(where(is.numeric)))) %>%
  dplyr::select(genes, total_counts, everything()) %>%
  dplyr::arrange(genes) %>%
  as.data.frame() -> d1

c2a %>%
  dplyr::mutate(total_counts = rowSums(across(where(is.numeric)))) %>%
  dplyr::select(genes, total_counts, everything()) %>%
  dplyr::arrange(genes) %>%
  as.data.frame() -> d2a
c3a %>%
  dplyr::mutate(total_counts = rowSums(across(where(is.numeric)))) %>%
  dplyr::select(genes, total_counts, everything()) %>%
  dplyr::arrange(genes) %>%
  as.data.frame() -> d3a
c4a %>%
  dplyr::mutate(total_counts = rowSums(across(where(is.numeric)))) %>%
  dplyr::select(genes, total_counts, everything()) %>%
  dplyr::arrange(genes) %>%
  as.data.frame() -> d4a

c2b %>%
  dplyr::mutate(total_counts = rowSums(across(where(is.numeric)))) %>%
  dplyr::select(genes, total_counts, everything()) %>%
  dplyr::arrange(genes) %>%
  as.data.frame() -> d2b
c3b %>%
  dplyr::mutate(total_counts = rowSums(across(where(is.numeric)))) %>%
  dplyr::select(genes, total_counts, everything()) %>%
  dplyr::arrange(genes) %>%
  as.data.frame() -> d3b
c4b %>%
  dplyr::mutate(total_counts = rowSums(across(where(is.numeric)))) %>%
  dplyr::select(genes, total_counts, everything()) %>%
  dplyr::arrange(genes) %>%
  as.data.frame() -> d4b

## build excel tables
cl1 <- list(d1, d2a, d3a, d4a, d2b, d3b, d4b)
names(cl1) <- c("SC_DAPI_counts",
                "SN1_DAPI_counts", "SN1_ChAT_counts", "SN1_NFH_counts",
                "SN2_DAPI_counts", "SN2_ChAT_counts", "SN2_NFH_counts")
writexl::write_xlsx(x = cl1, path = paste0(
  fp0, "new_filtered_SN_SC_DAPI_ChAT_NFH_pseudobulk_cts.xlsx"))

## save count matrices
saveRDS(object = d1, file = paste0(fp0, "SC_DAPI_counts.rds"))
saveRDS(object = d2a, file = paste0(fp0, "SN1_DAPI_counts.rds"))
saveRDS(object = d3a, file = paste0(fp0, "SN1_ChAT_counts.rds"))
saveRDS(object = d4a, file = paste0(fp0, "SN1_NFH_counts.rds"))
saveRDS(object = d2b, file = paste0(fp0, "SN2_DAPI_counts.rds"))
saveRDS(object = d3b, file = paste0(fp0, "SN2_ChAT_counts.rds"))
saveRDS(object = d4b, file = paste0(fp0, "SN2_NFH_counts.rds"))

################################# "gene plots" #################################
d2a %>%
  dplyr::select(-total_counts) %>%
  dplyr::full_join((d3a %>% dplyr::select(-total_counts)), by = "genes") %>%
  dplyr::full_join((d4a %>% dplyr::select(-total_counts)), by = "genes") %>%
  dplyr::full_join((d2b %>% dplyr::select(-total_counts)), by = "genes") %>%
  dplyr::full_join((d3b %>% dplyr::select(-total_counts)), by = "genes") %>%
  dplyr::full_join((d4b %>% dplyr::select(-total_counts)), by = "genes") %>%
  tidyr::gather(key = "sample", value = "counts", -genes) %>%
  dplyr::mutate(control = dplyr::case_when(
    stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
    stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
    stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
    stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43")) %>%
  dplyr::mutate(batch = dplyr::case_when(
    stringr::str_detect(sample, pattern = "SN1") ~ "SN1",
    stringr::str_detect(sample, pattern = "SN2") ~ "SN2")) %>%
  dplyr::mutate(stain = dplyr::case_when(
    stringr::str_detect(sample, pattern = "DAPI") ~ "DAPI",
    stringr::str_detect(sample, pattern = "ChAT") ~ "ChAT",
    stringr::str_detect(sample, pattern = "NFH") ~ "NFH")) %>%
  dplyr::mutate(control = as.factor(control),
                control = forcats::fct_relevel(control,
                                               c("CtrlFUSH+CtrlTDP43",
                                                 "FUSH", "TDP43",
                                                 "CtrlSOD1", "SOD1"))) %>%
  dplyr::arrange(genes) %>%
  as.data.frame() -> raw_counts


d1 %>%
  dplyr::select(-total_counts) %>%
  tidyr::gather(key = "sample", value = "counts", -genes) %>%
  dplyr::mutate(control = dplyr::case_when(
    stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
    stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
    stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
    stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43")) %>%
  dplyr::mutate(stain = dplyr::case_when(
    stringr::str_detect(sample, pattern = "DAPI") ~ "DAPI")) %>%
  dplyr::mutate(control = as.factor(control),
                control = forcats::fct_relevel(control,
                                               c("CtrlFUSH+CtrlTDP43",
                                                 "FUSH", "TDP43",
                                                 "CtrlSOD1", "SOD1"))) %>%
  dplyr::arrange(genes) %>%
  as.data.frame() -> raw_counts2

## plot the aggregated amounts
fgl1 <- list()
for (i in 1:length(unique(raw_counts$genes))) {
  raw_counts %>%
    dplyr::filter(stain == "DAPI",
                  genes == unique(raw_counts$genes)[i]) %>%
    ggplot(aes(x = control, y = counts, color = control)) +
    geom_point() +
    ggtitle(paste0(gsub("_", " ", unique(raw_counts$genes)[i]),
                   " non-normalized QC filtered SN DAPI pseudobulk counts")) +
    theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1),
          legend.position = "none") +
    facet_grid(~batch) -> fgl1[[i]]
}

fgl2 <- list()
for (i in 1:length(unique(raw_counts$genes))) {
  raw_counts %>%
    dplyr::filter(stain == "ChAT",
                  genes == unique(raw_counts$genes)[i]) %>%
    ggplot(aes(x = control, y = counts, color = control)) +
    geom_point() +
    ggtitle(paste0(gsub("_", " ", unique(raw_counts$genes)[i]),
                   " non-normalized QC filtered SN ChAT pseudobulk counts")) +
    theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1),
          legend.position = "none") +
    facet_grid(~batch) -> fgl2[[i]]
}

fgl3 <- list()
for (i in 1:length(unique(raw_counts$genes))) {
  raw_counts %>%
    dplyr::filter(stain == "NFH",
                  genes == unique(raw_counts$genes)[i]) %>%
    ggplot(aes(x = control, y = counts, color = control)) +
    geom_point() +
    ggtitle(paste0(gsub("_", " ", unique(raw_counts$genes)[i]),
                   " non-normalized QC filtered SN NFH pseudobulk counts")) +
    theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1),
          legend.position = "none") +
    facet_grid(~batch) -> fgl3[[i]]
}

fgl4 <- list()
for (i in 1:length(unique(raw_counts2$genes))) {
  raw_counts2 %>%
    dplyr::filter(stain == "DAPI",
                  genes == unique(raw_counts2$genes)[i]) %>%
    ggplot(aes(x = control, y = counts, color = control)) +
    geom_point() +
    ggtitle(paste0(gsub("_", " ", unique(raw_counts2$genes)[i]),
                   " non-normalized QC filtered SC DAPI pseudobulk counts")) +
    theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1),
          legend.position = "none") -> fgl4[[i]]
}

## save pseudobulk_gene_counts
pdf(file = paste0(fp0, "new_filtered_SN_DAPI_pseudobulk_cts.pdf"))
print(fgl1)
dev.off()

pdf(file = paste0(fp0, "new_filtered_SN_ChAT_pseudobulk_cts.pdf"))
print(fgl2)
dev.off()

pdf(file = paste0(fp0, "new_filtered_SN_NFH_pseudobulk_cts.pdf"))
print(fgl3)
dev.off()

pdf(file = paste0(fp0, "new_filtered_SC_DAPI_pseudobulk_cts.pdf"))
print(fgl4)
dev.off()

################################# PCAs #########################################
##### PCA SC DAPI
d1 %>%
  dplyr::arrange(desc(genes)) %>%
  dplyr::slice(-(1:2)) %>%
  dplyr::select(-total_counts) %>% 
  tibble::column_to_rownames(var = "genes") %>%
  t() %>% as.matrix() %>%
  prcomp(scale. = TRUE, center = TRUE) -> pca_1
pca_1$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sample") %>%
  dplyr::mutate(control = dplyr::case_when(
    stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
    stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
    stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
    stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43")) %>%
  dplyr::mutate(control = as.factor(control),
                control = forcats::fct_relevel(control,
                                               c("CtrlFUSH+CtrlTDP43",
                                                 "FUSH", "TDP43",
                                                 "CtrlSOD1", "SOD1"))) -> df1
autoplot(pca_1, data = df1, colour = 'control') +
  ggtitle("PCA of non-normalized QC filtered SC DAPI") -> p1

###### PCA SN1 DAPI
d2a %>%
  dplyr::arrange(desc(genes)) %>%
  dplyr::slice(-(1:2)) %>%
  dplyr::select(-total_counts) %>% 
  tibble::column_to_rownames(var = "genes") %>%
  t() %>% as.matrix() %>%
  prcomp(scale. = TRUE, center = TRUE) -> pca_2a
pca_2a$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sample") %>%
  dplyr::mutate(control = dplyr::case_when(
    stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
    stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
    stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
    stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43")) %>%
  dplyr::mutate(control = as.factor(control),
                control = forcats::fct_relevel(control,
                                               c("CtrlFUSH+CtrlTDP43",
                                                 "FUSH", "TDP43",
                                                 "CtrlSOD1", "SOD1"))) -> df2a
autoplot(pca_2a, data = df2a, colour = 'control') +
  ggtitle("PCA of non-normalized QC filtered SN1 DAPI") -> p2a

###### PCA SN1 Chat
d3a %>%
  dplyr::arrange(desc(genes)) %>%
  dplyr::slice(-(1:2)) %>%
  dplyr::select(-total_counts) %>% 
  tibble::column_to_rownames(var = "genes") %>%
  t() %>% as.matrix() %>%
  prcomp(scale. = TRUE, center = TRUE) -> pca_3a
pca_3a$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sample") %>%
  dplyr::mutate(control = dplyr::case_when(
    stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
    stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
    stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
    stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43")) %>%
  dplyr::mutate(control = as.factor(control),
                control = forcats::fct_relevel(control,
                                               c("CtrlFUSH+CtrlTDP43",
                                                 "FUSH", "TDP43",
                                                 "CtrlSOD1", "SOD1"))) -> df3a
autoplot(pca_3a, data = df3a, colour = 'control') +
  ggtitle("PCA of non-normalized QC filtered SN1 ChAT") -> p3a

###### PCA SN1 NFH
d4a %>%
  dplyr::arrange(desc(genes)) %>%
  dplyr::slice(-(1:2)) %>%
  dplyr::select(-total_counts) %>% 
  tibble::column_to_rownames(var = "genes") %>%
  t() %>% as.matrix() %>%
  prcomp(scale. = TRUE, center = TRUE) -> pca_4a
pca_4a$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sample") %>%
  dplyr::mutate(control = dplyr::case_when(
    stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
    stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
    stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
    stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43")) %>%
  dplyr::mutate(control = as.factor(control),
                control = forcats::fct_relevel(control,
                                               c("CtrlFUSH+CtrlTDP43",
                                                 "FUSH", "TDP43",
                                                 "CtrlSOD1", "SOD1"))) -> df4a
autoplot(pca_4a, data = df4a, colour = 'control') +
  ggtitle("PCA of non-normalized QC filtered SN1 NFH") -> p4a

###### PCA SN2 DAPI
d2b %>%
  dplyr::arrange(desc(genes)) %>%
  dplyr::slice(-(1:2)) %>%
  dplyr::select(-total_counts) %>% 
  tibble::column_to_rownames(var = "genes") %>%
  t() %>% as.matrix() %>%
  prcomp(scale. = TRUE, center = TRUE) -> pca_2b
pca_2b$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sample") %>%
  dplyr::mutate(control = dplyr::case_when(
    stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
    stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
    stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
    stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43")) %>%
  dplyr::mutate(control = as.factor(control),
                control = forcats::fct_relevel(control,
                                               c("CtrlFUSH+CtrlTDP43",
                                                 "FUSH", "TDP43",
                                                 "CtrlSOD1", "SOD1"))) -> df2b
autoplot(pca_2b, data = df2b, colour = 'control') +
  ggtitle("PCA of non-normalized QC filtered SN2 DAPI") -> p2b

###### PCA SN2 Chat #@ which(apply(h1 , 2, var)==0)
d3b %>%
  dplyr::arrange(desc(genes)) %>%
  dplyr::slice(-(1:2)) %>%
  dplyr::select(-total_counts) %>%
  dplyr::arrange(genes) %>%
  dplyr::filter(genes != "Apln") %>%
  tibble::column_to_rownames(var = "genes") %>%
  t() %>% as.matrix() %>%
  prcomp(scale. = TRUE, center = TRUE) -> pca_3b
pca_3b$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sample") %>%
  dplyr::mutate(control = dplyr::case_when(
    stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
    stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
    stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
    stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43")) %>%
  dplyr::mutate(control = as.factor(control),
                control = forcats::fct_relevel(control,
                                               c("CtrlFUSH+CtrlTDP43",
                                                 "FUSH", "TDP43",
                                                 "CtrlSOD1", "SOD1"))) -> df3b
autoplot(pca_3b, data = df3b, colour = 'control') +
  ggtitle("PCA of non-normalized QC filtered SN2 ChAT") -> p3b

###### PCA SN2 NFH
d4b %>%
  dplyr::arrange(desc(genes)) %>%
  dplyr::slice(-(1:2)) %>%
  dplyr::select(-total_counts) %>%
  dplyr::filter(genes != "Apln") %>%
  tibble::column_to_rownames(var = "genes") %>%
  t() %>% as.matrix() %>%
  prcomp(scale. = TRUE, center = TRUE) -> pca_4b
pca_4b$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sample") %>%
  dplyr::mutate(control = dplyr::case_when(
    stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
    stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
    stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
    stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43")) %>%
  dplyr::mutate(control = as.factor(control),
                control = forcats::fct_relevel(control,
                                               c("CtrlFUSH+CtrlTDP43",
                                                 "FUSH", "TDP43",
                                                 "CtrlSOD1", "SOD1"))) -> df4b
autoplot(pca_4b, data = df4b, colour = 'control') +
  ggtitle("PCA of non-normalized QC filtered SN2 NFH") -> p4b

############################ PCA SN1 + SN2 #####################################
###### PCA SN1 DAPI
d2a %>%
  dplyr::arrange(desc(genes)) %>%
  dplyr::slice(-(1:2)) %>%
  dplyr::select(-total_counts) %>%
  dplyr::left_join((d2b %>%
                       dplyr::arrange(desc(genes)) %>%
                       dplyr::slice(-(1:2)) %>%
                       dplyr::select(-total_counts)), by = "genes")  %>%
  tibble::column_to_rownames(var = "genes") %>%
  t() %>% as.matrix() %>%
  prcomp(scale. = TRUE, center = TRUE) -> pca_c2
pca_c2$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sample") %>%
  dplyr::mutate(control = dplyr::case_when(
    stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
    stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
    stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
    stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43")) %>%
  dplyr::mutate(control = as.factor(control),
                control = forcats::fct_relevel(control,
                                               c("CtrlFUSH+CtrlTDP43",
                                                 "FUSH", "TDP43",
                                                 "CtrlSOD1", "SOD1"))) %>%
  dplyr::mutate(batch = dplyr::case_when(
    stringr::str_detect(sample, pattern = "SN1") ~ "SN1",
    stringr::str_detect(sample, pattern = "SN2") ~ "SN2")) -> df2
autoplot(pca_c2, data = df2, colour = 'control', shape = "batch") +
  ggtitle("PCA of non-normalized QC filtered SN1 and SN2 DAPI") -> p2c

###### PCA SN1 Chat
d3a %>%
  dplyr::arrange(desc(genes)) %>%
  dplyr::slice(-(1:2)) %>%
  dplyr::select(-total_counts) %>%
  dplyr::left_join((d3b %>%
                       dplyr::arrange(desc(genes)) %>%
                       dplyr::slice(-(1:2)) %>%
                       dplyr::select(-total_counts)), by = "genes")  %>%
  tibble::column_to_rownames(var = "genes") %>%
  t() %>% as.matrix() %>%
  prcomp(scale. = TRUE, center = TRUE) -> pca_c3
pca_c3$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sample") %>%
  dplyr::mutate(control = dplyr::case_when(
    stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
    stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
    stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
    stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43")) %>%
  dplyr::mutate(control = as.factor(control),
                control = forcats::fct_relevel(control,
                                               c("CtrlFUSH+CtrlTDP43",
                                                 "FUSH", "TDP43",
                                                 "CtrlSOD1", "SOD1"))) %>%
  dplyr::mutate(batch = dplyr::case_when(
    stringr::str_detect(sample, pattern = "SN1") ~ "SN1",
    stringr::str_detect(sample, pattern = "SN2") ~ "SN2")) -> df3
autoplot(pca_c3, data = df3, colour = 'control', shape = "batch") +
  ggtitle("PCA of non-normalized QC filtered SN1 and SN2 Chat") -> p3c

###### PCA SN1 NFH
d4a %>%
  dplyr::arrange(desc(genes)) %>%
  dplyr::slice(-(1:2)) %>%
  dplyr::select(-total_counts) %>%
  dplyr::left_join((d4b %>%
                       dplyr::arrange(desc(genes)) %>%
                       dplyr::slice(-(1:2)) %>%
                       dplyr::select(-total_counts)), by = "genes")  %>%
  tibble::column_to_rownames(var = "genes") %>%
  t() %>% as.matrix() %>%
  prcomp(scale. = TRUE, center = TRUE) -> pca_c4
pca_c4$x %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "sample") %>%
  dplyr::mutate(control = dplyr::case_when(
    stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
    stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
    stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
    stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
    stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43")) %>%
  dplyr::mutate(control = as.factor(control),
                control = forcats::fct_relevel(control,
                                               c("CtrlFUSH+CtrlTDP43",
                                                 "FUSH", "TDP43",
                                                 "CtrlSOD1", "SOD1"))) %>%
  dplyr::mutate(batch = dplyr::case_when(
    stringr::str_detect(sample, pattern = "SN1") ~ "SN1",
    stringr::str_detect(sample, pattern = "SN2") ~ "SN2")) -> df4
autoplot(pca_c4, data = df4, colour = 'control', shape = "batch") +
  ggtitle("PCA of non-normalized QC filtered SN1 and SN2 NFH") -> p4c

pdf(file = paste0(fp0, "new_filtered_SN1_SN2_pseudobulk_PCAs.pdf"))
print(p1);
print(p2a); print(p3a); print(p4a); print(p2b); print(p3b); print(p4b)
print(p2c); print(p3c); print(p4c)
dev.off()

################################# "old code" ###################################

## plot gene loadings
# pca_res$x %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "genes") %>%
#   dplyr::arrange(desc(abs(PC2))) %>%
#   dplyr::mutate(genes = factor(genes, levels = unique(genes))) %>%
#   as.data.frame() -> pca_loadings

  # geom_text_repel(aes(label = gsub("_DAPI_SN1", "", sample))) +
  # theme(legend.position = "none") +

## make metadata table
# data.frame(file_name = as.factor(lf),
#            sample = paste0(gsub("-", "_", gsub("_labeled_count-matrix.csv", "", lf)), "_trans"),
#            stain = "DAPI",
#            slice = "trans",
#            stain_slice = "DAPI_trans",
#            group = as.factor(gsub("_.*", "", lf)),
#            control = as.factor(c(rep("CtrlFUSH+CtrlTDP43", 2), rep("CtrlSOD1", 3),
#                                  rep("CtrlFUSH+CtrlTDP43", 2), rep("FUSH", 3),
#                                  rep("SOD1", 3), rep("TDP43", 3)))) -> m2
# rownames(m2) <- m2$sample
# 
# ## run DESeq2 on "group"
# DESeqDataSetFromMatrix(countData = as.matrix(c2),
#                        colData = m2,
#                        design = ~ control) -> dds2
# keep2 <- rowSums(counts(dds2)) > 1
# dds2 <- DESeq(dds2[keep2, ])
# # resultsNames(dds2)
# rld2 <- rlogTransformation(dds2)
# PCA2 <- plotPCA(rld2, intgroup = "group") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   # geom_text_repel(aes(label = colnames(rld)),
#   #                 arrow = arrow(length = unit(0.03, "npc"),
#   #                 type = "closed", ends = "first"), force = 5) +
#   ggtitle(paste0("Principal component analysis of ", gsub("_", " ", proj2),
#                  " bulk counts"))

## make metadata table
# data.frame(file_name = as.factor(lf),
#            sample = paste0(gsub("-", "_", gsub(
#              "_labeled_non-expanded_count-matrix.csv", "", lf)), "_long"),
#            stain = "ChAT",
#            slice = "long",
#            stain_slice = "ChAT_long",
#            group = as.factor(gsub("_.*", "", lf)),
#            control = as.factor(c(rep("CtrlFUSH+CtrlTDP43", 2), rep("CtrlSOD1", 3),
#                                  rep("CtrlFUSH+CtrlTDP43", 2), rep("FUSH", 3),
#                                  rep("SOD1", 3), rep("TDP43", 3)))) -> m3
# rownames(m3) <- m3$sample
# 
# ## run DESeq2 on "group"
# DESeqDataSetFromMatrix(countData = as.matrix(c3),
#                        colData = m3,
#                        design = ~ control) -> dds3
# keep3 <- rowSums(counts(dds3)) > 1
# dds3 <- DESeq(dds3[keep3, ])
# # resultsNames(dds3)
# rld3 <- rlogTransformation(dds3)
# PCA3 <- plotPCA(rld3, intgroup = "group") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   # geom_text_repel(aes(label = colnames(rld)),
#   #                 arrow = arrow(length = unit(0.03, "npc"),
#   #                 type = "closed", ends = "first"), force = 5) +
#   ggtitle(paste0("Principal component analysis of ", gsub("_", " ", proj3),
#                  " bulk counts"))

## make metadata table
# data.frame(file_name = as.factor(lf),
#            sample = paste0(gsub("-", "_", gsub(
#              "_labeled_non-expanded_count-matrix.csv", "", lf)), "_long"),
#            stain = "NFH",
#            slice = "long",
#            stain_slice = "NFH_long",
#            group = as.factor(gsub("_.*", "", lf)),
#            control = as.factor(c(rep("CtrlFUSH+CtrlTDP43", 2), rep("CtrlSOD1", 3),
#                                  rep("CtrlFUSH+CtrlTDP43", 2), rep("FUSH", 3),
#                                  rep("SOD1", 3), rep("TDP43", 3)))) -> m4
# rownames(m4) <- m4$sample
# 
# ## run DESeq2 on "group"
# DESeqDataSetFromMatrix(countData = as.matrix(c4),
#                        colData = m4,
#                        design = ~ control) -> dds4
# keep4 <- rowSums(counts(dds4)) > 1
# dds4 <- DESeq(dds4[keep4, ])
# # resultsNames(dds4)
# rld4 <- rlogTransformation(dds4)
# PCA4 <- plotPCA(rld4, intgroup = "group") +
#   theme(plot.title = element_text(hjust = 0.5)) +
#   # geom_text_repel(aes(label = colnames(rld)),
#   #                 arrow = arrow(length = unit(0.03, "npc"),
#   #                 type = "closed", ends = "first"), force = 5) +
#   ggtitle(paste0("Principal component analysis of ", gsub("_", " ", proj4),
#                  " bulk counts"))

# h1 <- log10(c1)
# h1[is.infinite(h1)] <- 0

# pheatmap::pheatmap(mat = c1)

################################# "DAPI_long" ###########################
# proj1 <- "DAPI_longitudinal"
# fp <- "~/tmp/202305_piol/data/new_count_matrices/expanded/DAPI/long/"
# lf <- list.files(path = fp, pattern = "*.csv")
# dfs1 <- lapply(lf, function(i){readr::read_csv(file = paste0(fp, i))})
# dfs2 <- lapply(dfs1, function(i){i %>%
#     tibble::column_to_rownames(var = "cell_label")})
# dfs3 <- lapply(dfs2, function(i){t(i)})
# rl <- list()
# for (i in 1:length(dfs3)) {
#   rl[[i]] <- as.data.frame(rowSums(dfs3[[i]]))
#   print(nrow(rl[[i]]))
# }
# tibble::rownames_to_column(rl[[1]], var = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl[[2]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl[[3]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl[[4]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl[[5]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl[[6]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl[[7]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl[[8]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl[[9]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl[[10]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl[[11]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl[[12]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl[[13]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl[[14]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl[[15]], var = "genes"), by = "genes") %>%
#   dplyr::full_join(tibble::rownames_to_column(rl[[16]], var = "genes"), by = "genes") %>%
#   tibble::column_to_rownames(var = "genes") %>%
#   as.data.frame() -> c1
# c1[is.na(c1)] <- 0 ## coerce NA values to zero
# names(c1) <- c(paste0(gsub("_long", "", gsub("-", "_", gsub(
#   "_labeled_count-matrix.csv", "", lf))), "_long"))
# 
# ## make metadata table
# data.frame(file_name = as.factor(lf),
#            sample = paste0(gsub("_long", "", gsub("-", "_", gsub(
#              "_labeled_count-matrix.csv", "", lf))), "_long"),
#            stain = "DAPI", 
#            slice = "long",
#            stain_slice = "DAPI_long",
#            group = as.factor(gsub("_.*", "", lf)),
#            control = as.factor(c(rep("CtrlFUSH+CtrlTDP43", 2), rep("CtrlSOD1", 3),
#                                  rep("CtrlFUSH+CtrlTDP43", 2), rep("FUSH", 3),
#                                  rep("SOD1", 3), rep("TDP43", 3)))) -> m1
# rownames(m1) <- m1$sample

## run DESeq2 on "group"
# DESeqDataSetFromMatrix(countData = as.matrix(c1),
#                        colData = m1,
#                        design = ~ control) -> dds1
# keep1 <- rowSums(counts(dds1)) > 1
# dds1 <- DESeq(dds1[keep1, ])
# # resultsNames(dds1)
# rld1 <- rlogTransformation(dds1)
# PCA1 <- plotPCA(rld1, intgroup = "group") +
#         theme(plot.title = element_text(hjust = 0.5)) +
#         # geom_text_repel(aes(label = colnames(rld)),
#         #                 arrow = arrow(length = unit(0.03, "npc"),
#         #                 type = "closed", ends = "first"), force = 5) +
#         ggtitle(paste0("Principal component analysis of ", gsub("_", " ", proj1),
#                        " bulk counts"))


## reformat the count matrices
# c1 %>% 
#   tibble::rownames_to_column(var = "genes") %>%
#   dplyr::mutate(total_counts = rowSums(across(where(is.numeric)))) %>%
#   dplyr::select(genes, total_counts, everything()) %>%
#   dplyr::arrange(desc(total_counts)) %>%
#   as.data.frame() -> d1
# c2a %>% ##
#   # tibble::rownames_to_column(var = "genes") %>%
#   dplyr::mutate(total_counts = rowSums(across(where(is.numeric)))) %>%
#   dplyr::select(genes, total_counts, everything()) %>%
#   dplyr::arrange(desc(total_counts)) %>%
#   as.data.frame() -> d2a
# c3 %>%
#   tibble::rownames_to_column(var = "genes") %>%
#   dplyr::mutate(total_counts = rowSums(across(where(is.numeric)))) %>%
#   dplyr::select(genes, total_counts, everything()) %>%
#   dplyr::arrange(desc(total_counts)) %>%
#   as.data.frame() -> d3
# c4 %>%
#   tibble::rownames_to_column(var = "genes") %>%
#   dplyr::mutate(total_counts = rowSums(across(where(is.numeric)))) %>%
#   dplyr::select(genes, total_counts, everything()) %>%
#   dplyr::arrange(desc(total_counts)) %>%
#   as.data.frame() -> d4
# as.data.frame(counts(dds1, normalized = TRUE)) %>% 
#   tibble::rownames_to_column(var = "genes") %>%
#   dplyr::mutate(total_counts = rowSums(across(where(is.numeric)))) %>%
#   dplyr::select(genes, total_counts, everything()) %>%
#   dplyr::arrange(desc(total_counts)) %>%
#   as.data.frame() -> e1
# as.data.frame(counts(dds2, normalized = TRUE)) %>% 
#   tibble::rownames_to_column(var = "genes") %>%
#   dplyr::mutate(total_counts = rowSums(across(where(is.numeric)))) %>%
#   dplyr::select(genes, total_counts, everything()) %>%
#   dplyr::arrange(desc(total_counts)) %>%
#   as.data.frame() -> e2
# as.data.frame(counts(dds3, normalized = TRUE)) %>% 
#   tibble::rownames_to_column(var = "genes") %>%
#   dplyr::mutate(total_counts = rowSums(across(where(is.numeric)))) %>%
#   dplyr::select(genes, total_counts, everything()) %>%
#   dplyr::arrange(desc(total_counts)) %>%
#   as.data.frame() -> e3
# as.data.frame(counts(dds4, normalized = TRUE)) %>% 
#   tibble::rownames_to_column(var = "genes") %>%
#   dplyr::mutate(total_counts = rowSums(across(where(is.numeric)))) %>%
#   dplyr::select(genes, total_counts, everything()) %>%
#   dplyr::arrange(desc(total_counts)) %>%
#   as.data.frame() -> e4
# 
# ## build excel tables
# cl1 <- list(d1, e1,
#             d2, e2,
#             d3, e3,
#             d4, e4)
# names(cl1) <- c(gsub("_", " ", proj1), paste0(gsub("_", " ", proj1), " DESeq norm"),
#                 gsub("_", " ", proj2), paste0(gsub("_", " ", proj2), " DESeq norm"),
#                 gsub("_", " ", proj3), paste0(gsub("_", " ", proj3), " DESeq norm"),
#                 gsub("_", " ", proj4), paste0(gsub("_", " ", proj4), " DESeq norm"))
# 
# ## save count matrices
# saveRDS(object = c1, file = paste0(fp0, proj1, "_pseudobulk_count_matrix.rds"))
# saveRDS(object = c2, file = paste0(fp0, proj2, "_pseudobulk_count_matrix.rds"))
# saveRDS(object = c3, file = paste0(fp0, proj3, "_pseudobulk_count_matrix.rds"))
# saveRDS(object = c4, file = paste0(fp0, proj4, "_pseudobulk_count_matrix.rds"))
# writexl::write_xlsx(x = cl1, path = paste0(fp0, "DAPI_ChAT_NFH_pseudobulk_cts.xlsx"))
# 
# ## save PCAs
# pdf(file = paste0(fp0, "DAPI_ChAT_NFH_pseudobulk_PCAs.pdf"))
# print(PCA1); print(PCA2); print(PCA3); print(PCA4)
# dev.off()
# 
# ## genewise plots
# d1 %>% dplyr::select(-total_counts) %>% 
#   dplyr::full_join((d2 %>% dplyr::select(-total_counts)), by = "genes") %>% 
#   dplyr::full_join((d3 %>% dplyr::select(-total_counts)), by = "genes") %>% 
#   dplyr::full_join((d4 %>% dplyr::select(-total_counts)), by = "genes") %>%
#   dplyr::filter(!grepl("FP ", genes)) %>%
#   tidyr::gather(key = "sample", value = "counts", -genes) %>%
#   dplyr::mutate(control = dplyr::case_when(
#     stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
#     stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
#     stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
#     stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
#     stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
#     stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43")) %>%
#   dplyr::mutate(slice_type = dplyr::case_when(
#     stringr::str_detect(sample, pattern = "long") ~ "longitudinal",
#     stringr::str_detect(sample, pattern = "trans") ~ "transversal")) %>%
#   dplyr::mutate(stain_slice = dplyr::case_when(
#     stringr::str_detect(sample, pattern = "DAPI_long") ~ "DAPI_longitudinal",
#     stringr::str_detect(sample, pattern = "DAPI_trans") ~ "DAPI_transversal",
#     stringr::str_detect(sample, pattern = "ChAT") ~ "ChAT_transversal",
#     stringr::str_detect(sample, pattern = "NFH") ~ "NFH_transversal")) %>%
#   dplyr::mutate(stain_slice = as.factor(stain_slice),
#                 stain_slice = forcats::fct_relevel(stain_slice,
#                                         c("DAPI_longitudinal", "DAPI_transversal",
#                                           "ChAT_transversal", "NFH_transversal"))) %>%
#   dplyr::mutate(control = as.factor(control),
#                 control = forcats::fct_relevel(control,
#                                                c("CtrlFUSH+CtrlTDP43",
#                                                  "FUSH", "TDP43",
#                                                  "CtrlSOD1", "SOD1"))) %>%
#   as.data.frame() -> raw_counts
# 
# e1 %>% dplyr::select(-total_counts) %>% 
#   dplyr::full_join((e2 %>% dplyr::select(-total_counts)), by = "genes") %>% 
#   dplyr::full_join((e3 %>% dplyr::select(-total_counts)), by = "genes") %>% 
#   dplyr::full_join((e4 %>% dplyr::select(-total_counts)), by = "genes") %>%
#   dplyr::filter(!grepl("FP ", genes)) %>%
#   tidyr::gather(key = "sample", value = "counts", -genes) %>%
#   dplyr::mutate(control = dplyr::case_when(
#     stringr::str_detect(sample, pattern = "CtrlFUSH") ~ "CtrlFUSH+CtrlTDP43",
#     stringr::str_detect(sample, pattern = "CtrlSOD1") ~ "CtrlSOD1",
#     stringr::str_detect(sample, pattern = "CtrlTDP43") ~ "CtrlFUSH+CtrlTDP43",
#     stringr::str_detect(sample, pattern = "FUSH") ~ "FUSH",
#     stringr::str_detect(sample, pattern = "SOD1") ~ "SOD1",
#     stringr::str_detect(sample, pattern = "TDP43") ~ "TDP43")) %>%
#   dplyr::mutate(slice_type = dplyr::case_when(
#     stringr::str_detect(sample, pattern = "long") ~ "longitudinal",
#     stringr::str_detect(sample, pattern = "trans") ~ "transversal")) %>%
#   dplyr::mutate(stain_slice = dplyr::case_when(
#     stringr::str_detect(sample, pattern = "DAPI_long") ~ "DAPI_longitudinal",
#     stringr::str_detect(sample, pattern = "DAPI_trans") ~ "DAPI_transversal",
#     stringr::str_detect(sample, pattern = "ChAT") ~ "ChAT_transversal",
#     stringr::str_detect(sample, pattern = "NFH") ~ "NFH_transversal")) %>%
#   dplyr::mutate(stain_slice = as.factor(stain_slice),
#                 stain_slice = forcats::fct_relevel(stain_slice,
#                                                    c("DAPI_longitudinal", "DAPI_transversal",
#                                                      "ChAT_transversal", "NFH_transversal"))) %>%
#   dplyr::mutate(control = as.factor(control),
#                 control = forcats::fct_relevel(control,
#                                                c("CtrlFUSH+CtrlTDP43",
#                                                  "FUSH", "TDP43",
#                                                  "CtrlSOD1", "SOD1"))) %>%
#   as.data.frame() -> norm_counts
# 
# fgl1 <- list()
# for (i in 1:length(unique(raw_counts$genes))) {
#   raw_counts %>%
#     dplyr::filter(genes == unique(raw_counts$genes)[i]) %>%
#     ggplot(aes(x = control, y = counts, color = control)) +
#     geom_point() +
#     ggtitle(paste0(unique(raw_counts$genes)[i], " non-normalized pseudobulk counts")) +
#     theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1),
#           legend.position = "none") +
#     facet_grid(~stain_slice) -> fgl1[[i]]
# }
# 
# fgl2 <- list()
# for (i in 1:length(unique(norm_counts$genes))) {
#   raw_counts %>%
#     dplyr::filter(genes == unique(norm_counts$genes)[i]) %>%
#     ggplot(aes(x = control, y = counts, color = control)) +
#     geom_point() +
#     ggtitle(paste0(unique(raw_counts$genes)[i], " normalized pseudobulk counts")) +
#     theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust = 1),
#           legend.position = "none") +
#     facet_grid(~stain_slice) -> fgl2[[i]]
# }
# 
# ## save pseudobulk_gene_counts
# pdf(file = paste0(fp0, "DAPI_ChAT_NFH_pseudobulk_gene_counts.pdf"))
# print(fgl1)
# dev.off()
# pdf(file = paste0(fp0, "DAPI_ChAT_NFH_norm_pseudobulk_gene_counts.pdf"))
# print(fgl2)
# dev.off()

############################### old code #######################################
## run DESeq2 on "group"
# DESeqDataSetFromMatrix(countData = as.matrix(c1),
#                        colData = m1,
#                        design = ~ group) -> dds
# keep <- rowSums(counts(dds)) > 1
# dds <- DESeq(dds[keep, ])
# resultsNames(dds)
# rld <- rlogTransformation(dds)
# PCA <- plotPCA(rld, intgroup = "group") +
#         theme(plot.title = element_text(hjust = 0.5)) +
#         # geom_text_repel(aes(label = colnames(rld)),
#         #                 arrow = arrow(length = unit(0.03, "npc"),
#         #                 type = "closed", ends = "first"), force = 5) +
#         ggtitle("Principal component analysis")

## DE comparisons control FUS vs FUS
# res1 <- results(dds, contrast = c("group", "FUSH", "CtrlFUSH"), alpha = 0.05)
# merge(as.data.frame(res1), as.data.frame(counts(dds, normalized = TRUE)),
#       by = "row.names", sort = FALSE) %>%
#   as.data.frame() %>%
#   dplyr::rename(genes = names(.)[1]) %>%
#   dplyr::arrange(padj) -> rd1
# 
# ## DE comparisons control TDP vs TDP
# res2 <- results(dds, contrast = c("group", "TDP43", "CtrlTDP43"), alpha = 0.05)
# merge(as.data.frame(res2), as.data.frame(counts(dds, normalized = TRUE)),
#       by = "row.names", sort = FALSE) %>%
#   as.data.frame() %>%
#   dplyr::rename(genes = names(.)[1]) %>%
#   dplyr::arrange(padj) -> rd2
# 
# ## DE comparisons control SOD vs SOD
# res3 <- results(dds, contrast = c("group", "SOD1", "CtrlSOD1"), alpha = 0.05)
# merge(as.data.frame(res3), as.data.frame(counts(dds, normalized = TRUE)),
#       by = "row.names", sort = FALSE) %>%
#   as.data.frame() %>%
#   dplyr::rename(genes = names(.)[1]) %>%
#   dplyr::arrange(padj) -> rd3

## run DESeq2 on "control"
# DESeqDataSetFromMatrix(countData = as.matrix(c1),
#                        colData = m1,
#                        design = ~ control) -> dds2
# keep <- rowSums(counts(dds2)) > 1
# dds2 <- DESeq(dds2[keep, ])
# 
# ## DE comparisons control FUS & control TDP vs TDP
# res4 <- results(dds2, contrast = c("control", "TDP43", "Ctrl"), alpha = 0.05)
# merge(as.data.frame(res4), as.data.frame(counts(dds2, normalized = TRUE)),
#       by = "row.names", sort = FALSE) %>%
#   as.data.frame() %>%
#   dplyr::rename(genes = names(.)[1]) %>%
#   dplyr::arrange(padj) -> rd4
# 
# ## DE comparisons control FUS & control TDP vs FUS
# res5 <- results(dds2, contrast = c("control", "FUSH", "Ctrl2"), alpha = 0.05)
# merge(as.data.frame(res5), as.data.frame(counts(dds2, normalized = TRUE)),
#       by = "row.names", sort = FALSE) %>%
#   as.data.frame() %>%
#   dplyr::rename(genes = names(.)[1]) %>%
#   dplyr::arrange(padj) -> rd5
# 
# ## save results and wide count matrix
# c1 %>%  tibble::rownames_to_column(var = "genes") -> y1 
# res_list <- list(y1, rd1, rd2, rd3, rd4, rd5)
# names(res_list) <- c("raw pseudobulk counts", "ctrlFUSH vs FUSH",
#                      "ctrlTDP43 vs TDP43", "ctrlSOD1 vs SOD1",
#                      "ctrlFUSH and ctrlTDP43 vs TDP43",
#                      "ctrlFUSH and ctrlTDP43 vs FUSH")
# 
# # saveRDS(object = rd, file = paste0(fp, "TCGA_SKCM_DE_res.rds"))
# writexl::write_xlsx(x = res_list, path = paste0(fp, "piol_als_pseudobulk_res.xlsx"))
# 
# ## save PCA
# pdf(file = paste0(fp, "piol_als_pseudobulk_PCA.pdf"))
# print(PCA)
# dev.off()

# Can you also do the same but with these comparisons:
# - control FUS vs FUS
# - control TDP vs TDP
# - control FUS & control TDP vs TDP
# - control FUS & control TDP vs FUS
# - control SOD vs SOD