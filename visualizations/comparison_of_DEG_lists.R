## Gabriele's 10X Visium data analysis

# "DE_sc_Chat_pos_vs_neg_FUS.xlsx"  
# "DE_sc_Chat_pos_vs_neg_TDP43.xlsx"
# "DE_sc_Chat_pos_vs_neg.xlsx"

## I wanted to see the correlation of the DE genes in "DE_sc_Chat_pos_va_neg"
## with Nanostring "SC_chat_pos_vs_chat_neg" and all the datasets you did the
## correlation analysis and then do the same with the datasets TDP43 and FUS in
## the attachments and from Nanostring. So to do this last thing I need to send
## you the DE gene list for the comparison TDP_vs_control and FUS_vs_control of
## Nansotring data, but you can do the first part already.

## load libraries
library("dplyr")
library("tidyr")
library("readr")
library("tibble")
library("stringr")
library("ggplot2")
library("ggrepel")
library("DESeq2")
library("biomaRt")
library("readxl")
library("writexl")
library("janitor")
# library("DEGreport")
# library("RColorBrewer")
# library("pheatmap")

## file path
fp <- "~/Documents/tmp/2023/202308_piol/data/Gabriele_Visium/new_Gabriele_data/"
rp <- "~/Documents/tmp/2023/202308_piol/data/RDS/"
# lf <- list.files(fp, pattern = "*.xlsx")
proj <- "Gabriele_10X_Visium"
sig_qval = 0.05

## access Biomart API
# getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
#       mart = useDataset("mmusculus_gene_ensembl", useMart("ensembl"))) %>%
#   dplyr::rename(ensembl_id = ensembl_gene_id,
#                 genes = external_gene_name,
#                 entrez_id = entrezgene_id) %>% 
#   as.data.frame() -> mouse_biomart
# saveRDS(object = mouse_biomart, file = paste0(fp, "mouse_biomart.rds"))
mouse_biomart <- readRDS(file = paste0(fp, "mouse_biomart.rds"))
mouse_biomart %>%
  dplyr::filter(!duplicated(genes)) -> m1
dim(mouse_biomart); dim(m1)

## DE_sc_Chat_pos_vs_neg_FUS.xlsx DE results
readxl::read_excel(path = paste0(fp, "DE_sc_Chat_pos_vs_neg_FUS.xlsx"),
                   sheet = 1) %>% ## rank test
  as.data.frame() -> df1
readxl::read_excel(path = paste0(fp, "DE_sc_Chat_pos_vs_neg_FUS.xlsx"),
                   sheet = 2) %>% ## wald test
  dplyr::rename(genes = gene) %>%
  dplyr::left_join(m1, by = "genes") %>%
  dplyr::select(genes, ensembl_id, entrez_id, everything()) %>%
  dplyr::filter(qval < sig_qval) %>%
  as.data.frame() -> g2

## "DE_sc_Chat_pos_vs_neg_TDP43.xlsx" DE results
# readxl::read_excel(path = paste0(fp, "DE_sc_Chat_pos_vs_neg_TDP43.xlsx"),
#                    sheet = 1) %>%
#   as.data.frame() -> df3
readxl::read_excel(path = paste0(fp, "DE_sc_Chat_pos_vs_neg_TDP43.xlsx"),
                   sheet = 2) %>%
  dplyr::rename(genes = gene) %>%
  dplyr::left_join(m1, by = "genes") %>%
  dplyr::select(genes, ensembl_id, entrez_id, everything()) %>%
  dplyr::filter(qval < sig_qval) %>%
  as.data.frame() -> g4

## "DE_sc_Chat_pos_vs_neg.xlsx" DE results
# readxl::read_excel(path = paste0(fp, "DE_sc_Chat_pos_vs_neg.xlsx"),
#                    sheet = 1) %>%
#   as.data.frame() -> df5
readxl::read_excel(path = paste0(fp, "DE_sc_Chat_pos_vs_neg.xlsx"),
                   sheet = 2) %>%
  dplyr::rename(genes = gene) %>%
  dplyr::left_join(m1, by = "genes") %>%
  dplyr::select(genes, ensembl_id, entrez_id, everything()) %>%
  dplyr::filter(qval < sig_qval) %>%
  as.data.frame() -> g6

## collate list of lists
de_lists <- list(#df1,
                 g2,
                 #df3,
                 g4,
                 #df5,
                 g6)
names(de_lists) <- c(#"DE_sc_Chat_pos_vs_neg_FUS rank test",
                     "DE_sc_Chat_pos_vs_neg_FUS wald test",
                     #"DE_sc_Chat_pos_vs_neg_TDP43 rank test",
                     "DE_sc_Chat_pos_vs_neg_TDP43 wald test",
                     #"DE_sc_Chat_pos_vs_neg rank test",
                     "DE_sc_Chat_pos_vs_neg wald test")
# View(df1);
View(g2);
# View(df3);
View(g4);
# View(df5);
View(g6)

# saveRDS(object = g2, file = paste0(fp, "DE_sc_Chat_pos_vs_neg_FUS_wald_test.rds"))
# saveRDS(object = g4, file = paste0(fp, "DE_sc_Chat_pos_vs_neg_TDP43_wald_test.rds"))
# saveRDS(object = g6, file = paste0(fp, "DE_sc_Chat_pos_vs_neg_wald_test.rds"))

######################## comparing this data to other data sets ################
## rank by LFC
names(l0) <- c("Nanostring Chat+ vs Chat- SC",
               "Gitler GSE161621",
               "Levine_2 GSE190442",
               "Levine_1 GSE103892",
               "Le_Pichon GSE167597")
de_lists <- readRDS(file = paste0(rp, "all_DE_gene_lists.rds"))

## process nanostring
de_lists[[1]] %>%
  dplyr::select(-human_genes) %>%
  dplyr::rename(genes = Genes) %>%
  dplyr::left_join(mouse_biomart, by = "genes") %>%
  dplyr::select(genes, ensembl_id, entrez_id, everything()) -> ns_de

## process "Gitler GSE161621
de_lists[[2]] %>%
  dplyr::select(-c(group, contains("Gitler"))) %>%
  dplyr::left_join(mouse_biomart, by = "genes") %>%
  dplyr::select(genes, ensembl_id, entrez_id, everything()) -> gitler

## process "Levine_2 GSE190442"
readRDS(file = paste0("~/Documents/tmp/2023/202306_piol/DE_results/",
                                  "Levine_GSE190442_DE_RES2.rds")) %>%
  dplyr::select(-c(group, contains("Levine"))) %>%
  dplyr::rename(ensembl_id = ensembl_gene_id,
                synonym = external_synonym) %>%
  dplyr::left_join(mouse_biomart %>%
                     dplyr::select(-genes), by = "ensembl_id") %>%
  dplyr::select(human_genes, genes, ensembl_id, entrez_id, synonym, everything()
                ) -> levine_2

## process "Levine_1 GSE103892"
de_lists[[4]] %>%
  dplyr::select(-c(group, contains("Levine"))) %>%
  dplyr::left_join(mouse_biomart, by = "genes") %>%
  dplyr::select(genes, ensembl_id, entrez_id, everything()) -> levine_1

## process "Le_Pichon GSE167597"
de_lists[[5]] %>%
  dplyr::select(-c(group, contains("Pichon"))) %>%
  dplyr::left_join(mouse_biomart, by = "genes") %>%
  dplyr::select(genes, ensembl_id, entrez_id, everything()) -> le_pichon

## read formatted Gabriele's 10X Visium data analysis
g2 <- readRDS(file = paste0(fp, "DE_sc_Chat_pos_vs_neg_FUS_wald_test.rds"))
g4 <- readRDS(file = paste0(fp, "DE_sc_Chat_pos_vs_neg_TDP43_wald_test.rds"))
g6 <- readRDS(file = paste0(fp, "DE_sc_Chat_pos_vs_neg_wald_test.rds"))

## correct p-values
ns_de$pvalue[ns_de$pvalue == 0] <- 1e-300
ns_de$padj[ns_de$padj == 0] <- 1e-300
gitler$p_val[gitler$p_val == 0] <- 1e-300
gitler$p_val_adj[gitler$p_val_adj == 0] <- 1e-300
levine_2$p_val[levine_2$p_val == 0] <- 1e-300
levine_2$p_val_adj[levine_2$p_val_adj == 0] <- 1e-300
levine_1$p_val[levine_1$p_val == 0] <- 1e-300
levine_1$p_val[levine_1$p_val == 0] <- 1e-300
le_pichon$p_val[le_pichon$p_val == 0] <- 1e-300
le_pichon$p_val_adj[le_pichon$p_val_adj == 0] <- 1e-300

## build big DE list
big_de_list <- list(ns_de, gitler, levine_2, levine_1, le_pichon, g2, g4, g6)
names(big_de_list) <- c("Nanostring Chat+ vs Chat- SC",
                        "Gitler GSE161621", "Levine_2 GSE190442",
                        "Levine_1 GSE103892", "Le_Pichon GSE167597",
                        "DE_sc_Chat_pos_vs_neg_FUS",
                        "DE_sc_Chat_pos_vs_neg_TDP43",
                        "DE_sc_Chat_pos_vs_neg")
writexl::write_xlsx(x = big_de_list, path = paste0(rp, "big_DE_gene_lists.xlsx"))
saveRDS(object = big_de_list, file = paste0(rp, "big_DE_gene_lists.rds"))

## correlation
ns_de %>%
  dplyr::select(ensembl_id, log2FoldChange) %>%
  dplyr::rename(NS_LFC = log2FoldChange) %>%
  dplyr::filter(!duplicated(ensembl_id),
                !is.na(ensembl_id)) %>%
  dplyr::arrange(desc(NS_LFC)) %>%
  dplyr::left_join(g2 %>%
                     dplyr::select(ensembl_id, log2fc) %>%
                     dplyr::filter(!duplicated(ensembl_id),
                                   !is.na(ensembl_id)) %>%
                     dplyr::rename(DE_sc_Chat_pos_vs_neg_FUS_LFC = log2fc),
                   by = "ensembl_id") %>%
  dplyr::left_join(g4 %>%
                     dplyr::select(ensembl_id, log2fc) %>%
                     dplyr::filter(!duplicated(ensembl_id),
                                   !is.na(ensembl_id)) %>%
                     dplyr::rename(DE_sc_Chat_pos_vs_neg_TDP43_LFC = log2fc),
                   by = "ensembl_id") %>%
  dplyr::left_join(g6 %>%
                     dplyr::select(ensembl_id, log2fc) %>%
                     dplyr::filter(!duplicated(ensembl_id),
                                   !is.na(ensembl_id)) %>%
                     dplyr::rename(DE_sc_Chat_pos_vs_neg_LFC = log2fc),
                   by = "ensembl_id") %>%
  dplyr::left_join(levine_2 %>%
                     dplyr::select(ensembl_id, avg_log2FC) %>%
                     dplyr::filter(!duplicated(ensembl_id),
                                   !is.na(ensembl_id)) %>%
                     dplyr::rename(levine2_LFC = avg_log2FC),
                   by = "ensembl_id") %>%
  dplyr::left_join(le_pichon %>%
                     dplyr::select(ensembl_id, avg_log2FC) %>%
                     dplyr::filter(!duplicated(ensembl_id),
                                   !is.na(ensembl_id)) %>%
                     dplyr::rename(le_pichon_LFC = avg_log2FC),
                   by = "ensembl_id") %>%
  dplyr::left_join(gitler %>%
                     dplyr::select(ensembl_id, avg_log2FC) %>%
                     dplyr::filter(!duplicated(ensembl_id),
                                   !is.na(ensembl_id)) %>%
                     dplyr::rename(gitler_LFC = avg_log2FC),
                   by = "ensembl_id") %>%
  dplyr::left_join(levine_1 %>%
                     dplyr::select(ensembl_id, avg_log2FC) %>%
                     dplyr::filter(!duplicated(ensembl_id),
                                   !is.na(ensembl_id)) %>%
                     dplyr::rename(levine1_LFC = avg_log2FC),
                   by = "ensembl_id") %>%
  tibble::column_to_rownames(var = "ensembl_id") %>%
  as.matrix() -> z1

y1 <- rowSums(!is.na(z1)) > 2
x1 <- z1[y1,]
x1 %>%
  as.data.frame() %>%
  dplyr::select(-c(DE_sc_Chat_pos_vs_neg_FUS_LFC, DE_sc_Chat_pos_vs_neg_TDP43_LFC)) %>%
  dplyr::filter(abs(DE_sc_Chat_pos_vs_neg_LFC) < 200) %>%
  as.matrix() -> v1

writexl::write_xlsx(x = v1 %>%
                      as.data.frame() %>%
                      tibble::rownames_to_column(var = "ensembl_id"),
                    path = paste0(rp, "LFC_DE_gene_lists.xlsx"))

pheatmap::pheatmap(v1, scale = "none",
                   cluster_rows = FALSE,
                   cluster_cols = FALSE,
                   show_rownames = FALSE,
                   main = "LFC comparison across DE genes in motor neurons")

pheatmap::pheatmap(cor(v1, use = "na.or.complete"),
                   main = "Correlation of LFC across DE gene lists in motor neurons")

############################ Venn diagram time ################################
big_de_list <- readRDS(file = paste0(rp, "big_DE_gene_lists.rds"))