## diana instructions 

# Do these graphs:
# Volcano plot of the SC ChAT+ control samples vs SN ChAT+ control samples
#  
# And a new entry that I forgot before:
# Heatmap of these selected genes in SC ChAT+ control samples (4) vs SC ChAT-
# control samples (4):
data.frame(genes = c(
  "Pgk1", "Gapdh", "Hmbs", "Actb", "Chat", "Slc18a3", "Ahnak2", "Chodl",
  "Slc5a7", "Mmp9", "Pdgfd", "Grin3b", "Isl1", "Mnx1", "Vipr2", "Nefl", "Mapt",
  "Scn1a", "Kcna1", "Syn1", "Grin2a", "Grin1", "Syp", "Gria1", "Gria2", "Gfap",
  "Aldh1l1", "Aqp4", "Fgfr3", "S100b", "Mbp", "Sox10", "Mog", "Mag", "Cx3cr1",
  "Aif1", "Itgam", "Itgax", "Ptprcap")) -> cust_gl

## load libraries
library("dplyr")
library("tidyr")
library("tibble")
library("ggplot2")
library("ggrepel")
library("DESeq2")
# library("DEGreport")
library("RColorBrewer")
library("pheatmap")
library("readxl")
library("readr")
library("writexl")
library("janitor")

## file path
fp <- "~/Documents/tmp/dacruz/202307_piol/data/"
proj <- "SN_vs_SC"

############################## load input files ################################

## previous DE results
# readxl::read_excel(path = paste0(fp, "Results_SN_Vs_SC_Ordered.xlsx")) %>%
#   dplyr::select(Genes, 1:6) %>%
#   dplyr::arrange(padj) %>%
#   as.data.frame() -> de_res_snsc
# # View(de_res_old)
# 
# ## get top and bottom LFC genes
# dplyr::bind_rows(de_res_snsc %>%
#                   dplyr::arrange(log2FoldChange) %>%
#                   dplyr::filter(baseMean > 7) %>%
#                   dplyr::slice(1:20) %>%
#                   dplyr::arrange(desc(log2FoldChange)),
#                  de_res_snsc %>%
#                   dplyr::arrange(desc(log2FoldChange)) %>%
#                   dplyr::filter(baseMean > 7) %>%
#                   dplyr::slice(1:20)) -> top_LFC_genes
# 
readxl::read_excel(path = paste0(
  fp, "Maria_Count_Matrices/Normalized/Whole_Dataset/",
  "Count_Matrix_Normalized_Long.xlsx")) %>%
  as.data.frame() -> count_mat
# # View(count_mat)
# 
# count_mat %>%
#   dplyr::filter(!duplicated(sample)) -> meta_data
# 
# count_mat %>%
#   dplyr::filter(Genes %in% top_LFC_genes$Genes,
#                 class == "CTRL",
#                 segment == "ChATpos") %>%
#   dplyr::mutate(log_cts = log2(cts),
#                 Genes = factor(Genes, levels = unique(Genes)),
#                 Genes = forcats::fct_relevel(Genes, top_LFC_genes$Genes)
#                 ) -> snsc_sig
# is.na(snsc_sig) <- sapply(snsc_sig, is.infinite)
# snsc_sig[is.na(snsc_sig)] <- 0
# snsc_sig %>%
#   ggplot(aes(x = sample, y = Genes, fill = log_cts)) +
#   geom_tile() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
#   scale_fill_gradientn(colors = c("darkblue", "yellow", "red")) +
#   ggtitle(paste0("SN Chat+ vs SC Chat+ genes with highest DE LFC")) -> p1
# 
# snsc_sig %>%
#   dplyr::group_by(Genes) %>%
#   dplyr::mutate(zscore = (log_cts - mean(log_cts)) / sd(log_cts)) %>%
#   ggplot(aes(x = sample, y = Genes, fill = zscore)) +
#   geom_tile() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90)) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
#   ggtitle(paste0("SN Chat+ vs SC Chat+ genes with highest DE LFC by Z-score")
#           ) -> p2

## SN_Chatpos_vs_Chatneg actually significant
readRDS(file = paste0("~/Documents/tmp/dacruz/202307_piol/data/",
                      "SN_ChATpos_Vs_ChATneg_CTRL_metadata.rds")) %>% 
  dplyr::filter(class == "CTRL") %>%
  dplyr::arrange(desc(segment)) %>%
  as.data.frame() -> sn_meta

readxl::read_excel(path = paste0("~/Documents/tmp/dacruz/202306_piol/data/",
                                 "SN_Chatpos_vs_Chatneg_DEG_outliers.xlsx"),
                   sheet = 5) %>%
  dplyr::arrange(desc(padj)) %>%
  dplyr::select(genes, sn_meta$sample) %>%
  tidyr::gather(key = "sample", value = "cts", -genes) %>%
  dplyr::mutate(log_cts = log2(cts),
                sample = factor(sample, levels = unique(sample)),
                genes = factor(genes, levels = unique(genes))) %>%
  as.data.frame() -> sn_sig
is.na(sn_sig) <- sapply(sn_sig, is.infinite)
sn_sig[is.na(sn_sig)] <- 0
# View(sn_sig)

# sn_sig %>%
#   ggplot(aes(x = sample, y = genes, fill = log_cts)) +
#   geom_tile() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
#   scale_fill_gradientn(colors = c("darkblue", "yellow", "red")) +
#   ggtitle(paste0("SN Chat+ vs Chat- genes filtered significant DE LFC")) -> p3
# 
# sn_sig %>%
#   dplyr::group_by(genes) %>%
#   dplyr::mutate(zscore = (log_cts - mean(log_cts)) / sd(log_cts)) %>%
#   ggplot(aes(x = sample, y = genes, fill = zscore)) +
#   geom_tile() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
#   ggtitle(paste0(
#     "SN Chat+ vs Chat- genes filtered significant DE LFC by Z-score")) -> p4

# sn_sig %>% 
#   dplyr::left_join(sn_meta, by = "sample")  %>% 
#   dplyr::filter(genes == "Mpz") %>% 
#   ggplot(aes(x = sample, y = cts, color = segment)) +
#   geom_point()

## custom
count_mat %>%
  dplyr::filter(Genes %in% cust_gl$genes,
                class == "CTRL",
                region == "SC") %>%
  dplyr::mutate(log_cts = log2(cts),
                Genes = factor(Genes, levels = unique(Genes)),
                Genes = forcats::fct_relevel(Genes, rev(cust_gl$genes)))  %>%
  dplyr::arrange(desc(segment)) %>%
  dplyr::mutate(sample = factor(sample, levels = unique(sample)))  %>%
  dplyr::rename(genes = Genes) %>%
  as.data.frame() -> sn_sig3
is.na(sn_sig3) <- sapply(sn_sig3, is.infinite)
sn_sig3[is.na(sn_sig3)] <- 0
# sn_sig3 %>%
#   ggplot(aes(x = sample, y = genes, fill = log_cts)) +
#   geom_tile() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
#   scale_fill_gradientn(colors = c("darkblue", "yellow", "red")) +
#   ggtitle(paste0("SC Chat+ vs Chat- custom genes set")) -> p5

# sn_sig3 %>%
#   dplyr::group_by(genes) %>%
#   dplyr::mutate(zscore = (log_cts - mean(log_cts)) / sd(log_cts)) %>%
#   ggplot(aes(x = sample, y = genes, fill = zscore)) +
#   geom_tile() +
#   theme_classic() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 0)) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
#   ggtitle(paste0("SC Chat+ vs Chat- custom genes set by Z-score")) -> p6

# pdf(file = paste0(fp, "08082023_heatmaps.pdf"))
# p1; p2; p3; p4; p5; p6
# dev.off()

readxl::read_excel(path = paste0("~/Documents/tmp/dacruz/202306_piol/data/",
                                 "SN_Chatpos_vs_Chatneg_DEG_outliers.xlsx"),
                   sheet = 1) %>%
  dplyr::filter(#baseMean > 8.33,
                baseMean > 8,
                cooks_stat < 0.5,
                log2FoldChange > 0) %>%
  dplyr::arrange(desc(baseMean)) %>%
# dplyr::arrange(desc(padj)) %>%
# dplyr::select(genes, sn_meta$sample) %>%
# tidyr::gather(key = "sample", value = "cts", -genes) %>%
# dplyr::mutate(log_cts = log2(cts),
#               sample = factor(sample, levels = unique(sample)),
#               genes = factor(genes, levels = unique(genes))) %>%
  as.data.frame() -> sn_all

writexl::write_xlsx(x = sn_all, path = paste0("~/Documents/tmp/dacruz/202308_piol/",
                    "data/SN_Nanostring_reliably_detected_Chatpos_cutoff8.xlsx"))

# sn_all %>%
#   dplyr::filter(baseMean > 8,
#                 cooks_stat < 0.5,
#                 log2FoldChange > 0) %>%
#   dplyr::arrange(desc(baseMean)) %>%
#   dplyr::select(genes, sn_meta$sample[1:4]) %>%
#   tibble::column_to_rownames(var = "genes") %>%
#   as.matrix() -> df1
# mean(df1)


# writexl::write_xlsx(x = sn_all, path = paste0("~/Documents/tmp/dacruz/202308_piol/",
#                     "data/SN_Nanostring_reliably_detected.xlsx"))

## Volcanoes
# de_res_snsc %>% 
#   dplyr::filter(!is.na(padj)) %>% 
#   dplyr::mutate(threshold = as.factor(abs(log2FoldChange) > 1 & padj < 0.05)) %>%
#   as.data.frame() -> de_res_snsc2
# de_res_snsc2$padj[de_res_snsc2$padj == 0] <- 1e-300
# 
# de_res_snsc2 %>%
#   ggplot(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
#   geom_point(alpha = 0.4, size = 0.5) +
#   geom_hline(yintercept = -log10(0.05), color = "blue", linetype = "dashed") +
#   geom_vline(xintercept = -1, color = "blue", linetype = "dashed") +
#   geom_vline(xintercept = 1, color = "blue", linetype = "dashed") +
#   geom_text(label = de_res_snsc2$Genes, check_overlap = TRUE) +
#   theme(legend.position = "none") +
#   xlab("log2 fold change") + ylab("-log10 padj-value") +
#   scale_color_manual(values = c("#000000", "#FF0000")) +
#   ggtitle("Volcano plot of SN Chat+ vs SC Chat+ DE genes") -> p7
# 
# de_res_snsc2 %>%
#   ggplot(aes(x = log2FoldChange, y = -log10(padj), color = threshold)) +
#   geom_point(alpha = 0.4, size = 0.5) +
#   geom_hline(yintercept = -log10(0.05), color = "blue", linetype = "dashed") +
#   geom_vline(xintercept = -1, color = "blue", linetype = "dashed") +
#   geom_vline(xintercept = 1, color = "blue", linetype = "dashed") +
#   theme(legend.position = "none") +
#   xlab("log2 fold change") + ylab("-log10 padj-value") +
#   scale_color_manual(values = c("#000000", "#FF0000")) +
#   ggtitle("Volcano plot of SN Chat+ vs SC Chat+ DE genes") -> p8
# 
# pdf(file = paste0(fp, "SN+Chat+_vs_SC_Chat+_volcano_plots.pdf"))
# p7; p8
# dev.off()

######################################### old code #############################
# my_sample_col <- data.frame(segment =  as.character(meta_data$segment[meta_data$class == "CTRL"]))
# row.names(my_sample_col) <- colnames(a1)

# pheatmap::pheatmap(as.matrix(snsc_sig),
#                    # annotation_col = my_sample_col,
#                    main = paste("Heatmap of new significant DE genes for \n SN Chat+ vs Chat-"))
