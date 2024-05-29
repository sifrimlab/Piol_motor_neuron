### Make Lists of Reliable Genes with Volcano plots

# Diana:
# Starting from filtered data based on cook's distance and basemean >8.33 (1531
# genes that we called "reliably detected") generate volcano plots for the
# comparisons Chat-positive vs Chat-negative in control samples, and control vs
# FUSH in Chat-positive samples only. You can send me both volcano plots and
# excel files

## some references reflecting Maria DESeq2 parameters
## https://support.bioconductor.org/p/9135384/
## https://support.bioconductor.org/p/p132527/

## load libraries
library("dplyr")
library("tidyr")
library("tibble")
library("ggplot2")
library("DESeq2")
library("readr")
library("writexl")
library("pheatmap")
library("readxl")

## file paths
fp <- "~/Documents/tmp/dacruz/202307_piol/data/"
out <- "~/Documents/tmp/dacruz/202310_piol/data/"

## read long format count file
readxl::read_excel(path = paste0(
  fp, "Maria_Count_Matrices/Raw/Whole_Dataset/Count_Matrix_Long_Whole_Dataset.xlsx")) %>% 
  as.data.frame() -> raw_df

## convert long format dataframe to wide count matrix
raw_df %>% 
  dplyr::select(Genes, sample, cts) %>% 
  tidyr::pivot_wider(names_from = "Genes", values_from = "cts") %>%
  tibble::column_to_rownames(var = "sample") %>%
  t() %>% as.data.frame() -> count_mat

################# Chat-positive vs Chat-negative in control samples ############

## convert long data to metadata table (colData)
raw_df %>% 
  dplyr::filter(!duplicated(sample)) %>%
  dplyr::select(-c(Genes, cts, expected_neg)) %>%
  dplyr::mutate(group = paste0(region, "_", segment)) %>%
  dplyr::filter(class == "CTRL",
                region == "SN") %>%
  dplyr::mutate_all(as.factor) %>%
  as.data.frame() -> meta1
rownames(meta1) <- meta1$sample

## run DESeq2
dds <- DESeqDataSetFromMatrix(countData = count_mat %>%
                                dplyr::select(meta1$sample) %>%
                                as.matrix(),
                              colData = meta1,
                              design = ~ segment)
keep <- rowSums(counts(dds)) > 0
dds <- DESeq(dds[keep, ])
res <- results(dds, contrast = c("segment", "ChATpos", "ChATneg"), alpha = 0.05)

## this code chunk extracts the Cook's Distance
as.data.frame(mcols(dds)$maxCooks) %>%
  dplyr::rename(cooks_stat = names(.)[1]) %>%
  dplyr::arrange(desc(cooks_stat)) %>% 
  tibble::rownames_to_column(var = "genes") -> top_cooks

merge(as.data.frame(res), as.data.frame(counts(dds, normalized = TRUE)),
      by = "row.names", sort = FALSE) %>%
  as.data.frame() %>%
  dplyr::arrange(pvalue) %>%
  dplyr::rename(genes = names(.)[1]) %>%
  dplyr::left_join(top_cooks, by = "genes") %>% ## and this
  dplyr::select(genes, cooks_stat, everything()) -> res_data

## what is the average baseMean?
mean(res_data$baseMean)

## make "reliably detected" gene list
res_data %>% 
  dplyr::filter(baseMean > mean(res_data$baseMean),
                cooks_stat < 0.5) %>%
  dplyr::arrange(desc(baseMean)) -> rel_de

## how long is the "reliably detected" gene list
nrow(rel_de)

## save both lists as Excel
de_list1 <- list(res_data, rel_de)
names(de_list1) <- c("SN Chatpos vs Chatneg CTRL DEG", "reliably detected")
writexl::write_xlsx(x = de_list1, path = paste0(out, "SN_chatpos_vs_chatneg_ctrl_DEG.xlsx"))
# saveRDS(object = de_list1, file = paste0(out, "SN_chatpos_vs_chatneg_ctrl_DEG.rds"))

## volcano
pdf(file = paste0(out, "SN_chatpos_vs_chatneg_ctrl_DEG.pdf"))
rel_de %>%
  # dplyr::filter(!is.na(padj)) %>%
  dplyr::mutate(threshold = as.factor(abs(log2FoldChange) > 1 & pvalue < 0.05),
                group = as.factor(dplyr::case_when(
                  log2FoldChange > 1 & pvalue < 0.05 ~ "ChATpos",
                  log2FoldChange < -1 & pvalue < 0.05 ~ "ChATneg",
                  TRUE ~ "not significant"))) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue), color = group)) +
  geom_point(alpha = 0.5, size = 0.5) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_vline(xintercept = -1, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
  xlab("log2 fold change") +
  ylab("-log10 p-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("magenta", "limegreen", "gray3"))
dev.off()

##################### control vs FUSH in Chat-positive samples #################

## convert long data to metadata table (colData)
raw_df %>%
  dplyr::filter(!duplicated(sample)) %>%
  dplyr::select(-c(Genes, cts, expected_neg)) %>%
  dplyr::mutate(group = paste0(region, "_", segment)) %>%
  dplyr::filter(group == "SN_ChATpos",
                region == "SN",
                class != "TDP_43") %>%
  dplyr::mutate_all(as.factor) %>%
  as.data.frame() -> meta2
rownames(meta2) <- meta2$sample

## run DESeq2
dds <- DESeqDataSetFromMatrix(countData = count_mat %>%
                                dplyr::select(meta2$sample) %>%
                                as.matrix(),
                              colData = meta2,
                              design = ~ class)
keep <- rowSums(counts(dds)) > 0
dds <- DESeq(dds[keep, ])
res <- results(dds, contrast = c("class", "CTRL", "FUS"), alpha = 0.05)

barplot(dds$sizeFactor, main = "Size Factors for samples", col = "lavender",
        cex.axis = 0.7, las = 2, ylab = "Size Factors", cex.names = 0.4)

as.data.frame(mcols(dds)$maxCooks) %>%
  dplyr::rename(cooks_stat = names(.)[1]) %>%
  dplyr::arrange(desc(cooks_stat)) %>% 
  tibble::rownames_to_column(var = "genes") -> top_cooks

merge(as.data.frame(res), as.data.frame(counts(dds, normalized = TRUE)),
      by = "row.names", sort = FALSE) %>%
  as.data.frame() %>%
  dplyr::arrange(pvalue) %>%
  dplyr::rename(genes = names(.)[1]) %>%
  dplyr::left_join(top_cooks, by = "genes") %>%
  dplyr::select(genes, cooks_stat, everything()) -> res_data

## what is the average baseMean?
mean(res_data$baseMean)

## make "reliably detected" gene list
res_data %>% 
  dplyr::filter(baseMean > mean(res_data$baseMean),
                cooks_stat < 0.5) %>%
  dplyr::arrange(desc(baseMean)) -> rel_de

## how long is the "reliably detected" gene list
nrow(rel_de)

## save both lists as Excel
de_list2 <- list(res_data, rel_de)
names(de_list2) <- c("SN CTRL vs FUS Chatpos full DEG", "reliably detected")
writexl::write_xlsx(x = de_list2, path = paste0(out, "SN_CTRL_vs_FUS_in_chatpos_DEG.xlsx"))
# saveRDS(object = de_list2, file = paste0(out, "SN_CTRL_vs_FUS_in_chatpos_DEG.rds"))

## volcano
pdf(file = paste0(out, "SN_CTRL_vs_FUS_in_chatpos_DEG.pdf"))
rel_de %>%
  # dplyr::filter(!is.na(padj)) %>%
  dplyr::mutate(threshold = as.factor(abs(log2FoldChange) > 1 & pvalue < 0.05),
                group = as.factor(dplyr::case_when(
                  log2FoldChange > 1 & pvalue < 0.05 ~ "CTRL",
                  log2FoldChange < -1 & pvalue < 0.05 ~ "FUS",
                  TRUE ~ "not significant"))) %>%
  ggplot(aes(x = log2FoldChange, y = -log10(pvalue), color = group)) +
  geom_point(alpha = 0.75, size = 0.75) +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_vline(xintercept = -1, color = "red", linetype = "dashed") +
  geom_vline(xintercept = 1, color = "red", linetype = "dashed") +
  xlab("log2 fold change") +
  ylab("-log10 p-value") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = c("magenta", "limegreen", "gray3"))
dev.off()
