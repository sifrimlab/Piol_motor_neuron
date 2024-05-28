## this script compares the SN and SC regions and finds DE genes


## load libraries
library("dplyr")
library("tidyr")
library("tibble")
library("ggplot2")
library("ggrepel")
library("DESeq2")
library("DEGreport")
library("RColorBrewer")
library("pheatmap")
library("readxl")
library("readr")
library("writexl")
library("janitor")

## file path
fp <- "~/Documents/tmp/2023/202307_piol/data/"
proj <- "SN_SC_chat+"

############################## load input files ################################

## read count matrix
readxl::read_excel(path = paste0(
  fp, "Maria_Count_Matrices/Raw/Whole_Dataset/Count_Matrix_Long_Whole_Dataset.xlsx")) %>% 
  as.data.frame() -> raw_df

raw_df %>% 
  dplyr::select(Genes, sample, cts) %>% 
  tidyr::pivot_wider(names_from = "Genes", values_from = "cts") %>%
  tibble::column_to_rownames(var = "sample") %>%
  t() %>% as.data.frame() -> count_mat

## make metadata
raw_df %>% 
  dplyr::filter(!duplicated(sample)) %>% 
  dplyr::select(-c(Genes, cts, expected_neg)) %>%
  dplyr::mutate(group = paste0(region, "_", segment)) %>%
  dplyr::filter(class == "CTRL",
                segment == "ChATpos") %>%
  dplyr::mutate_all(as.factor) %>%
  as.data.frame() -> meta_data
rownames(meta_data) <- meta_data$sample

################################# run DESeq2 ###################################
dds <- DESeqDataSetFromMatrix(countData = count_mat %>%
                                            dplyr::select(meta_data$sample) %>%
                                            as.matrix(),
                              colData = meta_data,
                              design = ~ group)
keep <- rowSums(counts(dds)) > 0
dds <- DESeq(dds[keep, ])
res <- results(dds, contrast = c("group", "SN_ChATpos", "SC_ChATpos"),
               alpha = 0.05)
res %>% 
  as.data.frame() %>%
  tibble::rownames_to_column(var = "genes") %>%
  dplyr::arrange(padj) -> de_res
writexl::write_xlsx(x = de_res, path = paste0(fp, "SN_SC_chatpos_DE_res.xlsx"))
saveRDS(object = de_res, file = paste0(fp, "SN_SC_chatpos_DE_res.rds"))