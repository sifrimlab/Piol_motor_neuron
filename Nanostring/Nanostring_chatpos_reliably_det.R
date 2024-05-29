## reliably detected genes in SN Chat pos (only)
## mean counts divided by area

## load libraries
library("dplyr")
library("tidyr")
library("tibble")
library("ggplot2")
library("readxl")
library("readr")
library("writexl")
library("matrixStats")
# library("janitor")
# library("ggrepel")
# library("RColorBrewer")
# library("pheatmap")

## file path
fp <- "~/Documents/tmp/dacruz/202307_piol/data/"
proj <- "SN_chatpos"

## read nanostring metadata table
readxl::read_excel(path = paste0("~/Documents/tmp/dacruz/Documents/tmp/dacruz/202309_piol/data/",
                                 "Metadata_SN_DP.xlsx")) %>%
  dplyr::filter(class == "CTRL", ## only control samples
                segment == "ChATpos") %>%
  dplyr::mutate(area = as.double(area),
                area_mean = mean(area),
                area_ratio = area / area_mean) %>%
  as.data.frame() -> Metadata_SN_DP
# View(Metadata_SN_DP)

## read bulk normalized counts
readxl::read_excel(path = paste0("~/Documents/tmp/dacruz/202306_piol/data/",
                                 "SN_Chatpos_vs_Chatneg_DEG_outliers.xlsx"),
                   sheet = 1) %>%
  dplyr::select(genes, Metadata_SN_DP$sample) %>%
  dplyr::mutate(
    Chat_pos_ctrl1 = DSP_1005580000422_D_A02.dcc / Metadata_SN_DP$area_ratio[1],
    Chat_pos_ctrl2 = DSP_1005580000422_D_A10.dcc / Metadata_SN_DP$area_ratio[2],
    Chat_pos_ctrl3 = DSP_1005580000422_D_B10.dcc / Metadata_SN_DP$area_ratio[3],
    Chat_pos_ctrl4 = DSP_1005580000422_D_B12.dcc / Metadata_SN_DP$area_ratio[4]) %>%
  dplyr::select(genes, contains("Chat")) %>%
  tibble::column_to_rownames(var = "genes") %>%
  as.data.frame() -> sn_chat_pos_ctrl

## calculate the standard deviation for each row fo the normalized counts
transform(sn_chat_pos_ctrl, SD = rowSds(as.matrix(sn_chat_pos_ctrl), na.rm = TRUE)) %>%
  dplyr::mutate(mean_counts = rowMeans(across(1:4))) %>%
  dplyr::filter(mean_counts > 0) %>%
  tibble::rownames_to_column(var = "genes") %>%
  tidyr::gather(key = "samples", value = "counts", -c(genes, SD, mean_counts)) %>%
  dplyr::group_by(genes) %>%
  dplyr::mutate(zscore = (counts - mean(counts)) / sd(counts)) %>%
  as.data.frame() -> sn_chat_pos_ctrl_sd
# nrow(sn_chat_pos_ctrl_sd)

sn_chat_pos_ctrl_sd %>% 
  dplyr::arrange(desc(zscore)) %>%
  tibble::rowid_to_column(var = "ID") %>%
  ggplot(aes(zscore)) +
  geom_density() +
  # xlab("genes") +
  # geom_hline(yintercept = global_mean, color = "red", linetype = "dashed") +
  ggtitle("Distribution of mean gene counts in SN Chat+ samples")


## global mean count after area normalization
global_mean <- mean(as.matrix(sn_chat_pos_ctrl_sd[,1:4]), na.rm = TRUE)
global_mean

## plot mean gene count distribution
sn_chat_pos_ctrl_sd %>% 
  tibble::rownames_to_column(var = "gene") %>%
  dplyr::arrange(desc(mean_counts)) %>%
  tibble::rowid_to_column(var = "ID") %>%
  ggplot(aes(x = ID, y = mean_counts)) +
  geom_line() +
  xlab("genes") +
  geom_hline(yintercept = global_mean, color = "red", linetype = "dashed") +
  ggtitle("Distribution of mean gene counts in SN Chat+ samples")


# X %>% 
#   mutate(ID = row_number()) %>%
#   group_by(ID) %>%
#   do(data.frame(., SD = sd(unlist(.[vars_to_sum]), na.rm=T)))

# writexl::write_xlsx(x = sn_all,
#                     path = paste0("~/Documents/tmp/dacruz/202308_piol/data/",
#                                   "SN_Nanostring_reliably_detected_Chatpos_cutoff8.xlsx"))


# readxl::read_excel(path = paste0("~/Documents/tmp/dacruz/202306_piol/data/",
#                                  "SN_Chatpos_vs_Chatneg_DEG_outliers.xlsx"),
#                    sheet = 1) %>%
#   dplyr::filter(#baseMean > 8.33,
#     baseMean > 8,
#     cooks_stat < 0.5,
#     log2FoldChange > 0) %>%
#   dplyr::arrange(desc(baseMean)) %>%
#   # dplyr::arrange(desc(padj)) %>%
#   # dplyr::select(genes, sn_meta$sample) %>%
#   # tidyr::gather(key = "sample", value = "cts", -genes) %>%
#   # dplyr::mutate(log_cts = log2(cts),
#   #               sample = factor(sample, levels = unique(sample)),
#   #               genes = factor(genes, levels = unique(genes))) %>%
#   as.data.frame() -> sn_all