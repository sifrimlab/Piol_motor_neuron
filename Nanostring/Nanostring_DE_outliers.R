## outlier detection

## references:
## https://support.bioconductor.org/p/68277/
## https://rpubs.com/DragonflyStats/Cooks-Distance
## https://www.nonlinear.com/support/progenesis/comet/faq/v2.0/pq-values.aspx
## http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#pvaluesNA

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
proj <- "SN_outliers"

############################## load input files ################################

## read count matrix
readr::read_csv(file = paste0(
  fp, "Maria_Count_Matrices/Raw/SN/Count_Matrix_Raw_Long_SN.csv")) %>% 
  as.data.frame() -> raw_df

raw_df %>% 
  dplyr::select(Genes, sample, cts) %>% 
  tidyr::pivot_wider(names_from = "Genes", values_from = "cts") %>%
  tibble::column_to_rownames(var = "sample") %>%
  t() %>% as.data.frame() -> count_mat

## make metadata
raw_df %>% 
  dplyr::filter(!duplicated(sample)) %>% 
  dplyr::select(-c(Genes, cts, expected_neg, scan_name, region)) %>%
  dplyr::mutate(group = paste0(class, "_", segment)) %>%
  # dplyr::filter(class == "CTRL") %>%
  # dplyr::filter(sample != "DSP_1005580000422_D_B10.dcc") %>%
  dplyr::mutate_all(as.factor) %>%
  as.data.frame() -> meta_data
rownames(meta_data) <- meta_data$sample
# saveRDS(object = meta_data, file = paste0(fp, "SN_ChATpos_Vs_ChATneg_CTRL_metadata.rds"))

## previous DE results
readxl::read_excel(path = paste0(fp, "Results_SN_ChATpos_Vs_ChATneg_CTRL_Ordered.xlsx")) %>%
  dplyr::select(Genes, 1:6) %>%
  dplyr::arrange(pvalue) %>%
  as.data.frame() -> de_res_old
# View(de_res_old)

################################# run DESeq2 ###################################
dds <- DESeqDataSetFromMatrix(countData = count_mat %>% 
                                            dplyr::select(meta_data$sample) %>%
                                            as.matrix(),
                              colData = meta_data,
                              design = ~ group)
keep <- rowSums(counts(dds)) > 0
dds <- DESeq(dds[keep, ])
res <- results(dds, contrast = c("group", "CTRL_ChATpos", "CTRL_ChATneg"),
               alpha = 0.05)

## PCA
rld <- DESeq2::vst(dds)
plotPCA(rld, intgroup = "group") +
        theme(plot.title = element_text(hjust = 0.5)) +
        geom_text_repel(aes(label = colnames(rld)),
                        arrow = arrow(length = unit(0.03, "npc"),
                        type = "closed", ends = "first")) +
        ggtitle("PCA by genotype") -> p1

## process DE results
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

## sign DE genes 
res_data %>%
    dplyr::filter(pvalue < 0.05) -> res_data_sig

## 
res_data %>% 
  dplyr::filter(pvalue < 0.05,
                cooks_stat > 0.5) -> res_data_sig_outlier

## all outliers
res_data %>%
    dplyr::filter(cooks_stat > 0.5) -> res_data_all_outlier

res_data %>%
    dplyr::filter(padj < 0.05,
                  # abs(log2FoldChange) > 1,
                  cooks_stat < 0.5,
                  baseMean > 8.33
                  ) -> actually_sig

# list(res_data,
#      res_data_sig,
#      res_data_sig_outlier,
#      res_data_all_outlier,
#      actually_sig) -> de_list
# names(de_list) <- c("All DE genes",
#                     "Sig. DE genes",
#                     "Sig. DE genes - outliers",
#                     "Non-Sig. DE genes - outliers",
#                     "Actually Sig.")
# writexl::write_xlsx(x = de_list,
#                     path = paste0(fp, "SN_Chatpos_vs_Chatneg_DEG_outliers.xlsx"))
# 
# res_data %>% 
#   ggplot(aes(x = log10(baseMean))) +
#   geom_vline(xintercept = log10(8.33), color = "red", linetype = "dashed") +
#   geom_density() +
#   ggtitle("Density distribution of log10 baseMean normalized counts") -> p2
#   
# 
# actually_sig %>% 
#   dplyr::select(genes, contains(c("A02", "A03", "A10", "A11", "B10", "B11", "B12", "C01"))) %>% 
#   tibble::column_to_rownames(var = "genes") %>% 
#   log2() -> a1
# 
# is.na(a1) <- sapply(a1, is.infinite)
# a1[is.na(a1)] <- 0
# my_sample_col <- data.frame(segment =  as.character(meta_data$segment[meta_data$class == "CTRL"]))
# row.names(my_sample_col) <- colnames(a1)
# 
# pheatmap::pheatmap(as.matrix(a1),
#                    annotation_col = my_sample_col,
#                    main = paste("Heatmap of new significant DE genes for \n SN Chat+ vs Chat-"))

# summary(res)[1]
# 
# View(assays(dds)[["cooks"]])

################################# QC figures ###################################
# summary(res)

# plot(metadata(res)$filterNumRej, 
#      type="b", ylab="number of rejections",
#      xlab="quantiles of filter")
# lines(metadata(res)$lo.fit, col="red")
# abline(v=metadata(res)$filterTheta)


# par(mfrow = c(1, 1))
# boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2,
#         main = paste0("Boxplot of the Cook's distances for each Control sample"))
# abline(h = log10(4/nrow(count_mat)), col = "red")

# res_data %>% 
#   dplyr::filter(is.na(padj)) -> res_data_NA
# nrow(res_data_NA)

################################# get outliers  ################################

## get the cooks distance
# assays(dds)[["cooks"]] %>%
#   as.data.frame() %>%
#   tibble::rownames_to_column(var = "Genes") %>%
#   dplyr::arrange("Genes") -> cooks_dist1
# 
# assays(dds)[["cooks"]] %>%
#   as.data.frame() -> cooks_dist2
# 
# pheatmap::pheatmap(as.matrix(cooks_dist2))
# 
# apply(assays(dds)[["cooks"]], 1, max) %>%
#   as.data.frame() %>%
#   dplyr::rename(cook_stat = names(.)[1]) %>%
#   tibble::rownames_to_column(var = "Genes") -> max_cooks
# 
# max_cooks %>%
#   dplyr::filter(cook_stat < (4/(nrow(meta_data)))) -> pass_filter
# nrow(pass_filter) / nrow(res_data)


# 
# mcols(dds)$maxCooks <- apply(assays(dds)[["cooks"]], 1, max)
# plot(mcols(dds)$baseMean, mcols(dds)$maxCooks)
# # this requires you not filter or order the 'res' object
# stopifnot(all.equal(rownames(dds), rownames(res)))
# plot(res$log2FoldChange, mcols(dds)$maxCooks)
## After looking at examples and upon deciding you want to filter on Cook's
## distance per gene, you can specify cooksCutoff:

# res$pvalue[mcols(dds)$maxCooks > cooksCutoff] <- NA
# # optionally, also mean filtering:
# res$pvalue[res$baseMean < meanFilter] <- NA
# res$padj <- p.adjust(res$pvalue, method = "BH")

res_data %>%
  dplyr::arrange(log2FoldChange) %>%
  dplyr::slice(1:20) %>%
  dplyr::arrange(desc(log2FoldChange)) -> bottom_20

res_data %>%
  dplyr::arrange(desc(log2FoldChange)) %>%
  dplyr::slice(1:20) -> top_20