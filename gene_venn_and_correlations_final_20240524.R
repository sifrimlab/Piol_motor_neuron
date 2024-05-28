## establishing uniqueness in DE results by mapping to ensembl IDs

## load libraries
library("dplyr")
library("stringr")
library("tibble")
library("biomaRt")
library("grid")
library("writexl")
library("ggplot2")
library("corrplot")
library("ggpubr")
library("ggrepel")
library("gridExtra")
library("VennDiagram")

## Helper function to display Venn diagram
display_venn <- function(x, ...){
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}

## file path
fp <- "~/Documents/tmp/dacruz/202306_piol/"
rp <- "~/Documents/tmp/dacruz/202403_piol/data/"

## read Yadav/Levine 2 DE results
readRDS(file = paste0("~/Documents/tmp/dacruz/202306_piol/DE_results/RDS/GSE190442_DE.rds")) %>%
  dplyr::rename(genes = names(.)[1]) %>%
  dplyr::mutate(Levine_GSE190442_gene_rank = row_number()) %>%
  as.data.frame() -> Levine_GSE190442_DE_RES

## read Gautier/Gitler 2 DE results
readRDS("~/Documents/tmp/dacruz/202402_piol/Gauthier_data/GSE228778_motoneuron_DEGs.rds") %>%
  dplyr::rename(genes = names(.)[1]) %>%
  dplyr::mutate(Levine_GSE190442_gene_rank = row_number()) %>%
  as.data.frame() -> Gautier_GSE228778_DE_RES

## retrieve list of official mouse gene symbols
getBM(attributes = c("external_synonym", "ensembl_gene_id", "external_gene_name"),
      mart = useDataset("mmusculus_gene_ensembl",
                        useMart("ensembl", host = "https://www.ensembl.org"))) %>%
  dplyr::rename(genes = external_gene_name) %>%
  dplyr::filter(stringr::str_length(genes) > 1) %>% ## no blank gene entries
  as.data.frame() -> mouse_gene_ids

## get Marts and convert human gene list to mouse marker genes
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                host = "https://dec2021.archive.ensembl.org")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                host = "https://dec2021.archive.ensembl.org")
genesV2 = getLDS(attributes = "hgnc_symbol",
                 filters = "hgnc_symbol",
                 values = Levine_GSE190442_DE_RES$gene,
                 mart = human,
                 attributesL = "mgi_symbol",
                 martL = mouse,
                 uniqueRows = TRUE)

## map human to mouse genes
Levine_GSE190442_DE_RES %>% 
  dplyr::left_join((genesV2 %>%
                     dplyr::rename(genes = names(.)[1])),
                   by = "genes") %>% 
  dplyr::filter(!is.na(MGI.symbol),
                !duplicated(MGI.symbol)) %>% 
  dplyr::rename(human_genes = genes) %>%
  dplyr::rename(genes = MGI.symbol) %>%
  dplyr::left_join(mouse_gene_ids, by = "genes") %>%
  dplyr::filter(!duplicated(ensembl_gene_id),
                !is.na(ensembl_gene_id)) %>%
  as.data.frame() -> Levine_GSE190442_DE_RES2

## map human to mouse genes
Gautier_GSE228778_DE_RES %>%
  dplyr::left_join((genesV2 %>%
                      dplyr::rename(genes = names(.)[1])),
                   by = "genes") %>% 
  dplyr::filter(!is.na(MGI.symbol),
                !duplicated(MGI.symbol)) %>% 
  dplyr::rename(human_genes = genes) %>%
  dplyr::rename(genes = MGI.symbol) %>%
  dplyr::left_join(mouse_gene_ids, by = "genes") %>%
  dplyr::filter(!duplicated(ensembl_gene_id),
                !is.na(ensembl_gene_id)) %>%
  as.data.frame() -> Gautier_GSE228778_DE_RES2

## read Nanostring results
readRDS(file = paste0(fp, "DE_results/RDS/DE_Nanostring_sig_genes.rds")) %>%
  dplyr::rename(genes = names(.)[1]) %>%
  dplyr::select(-human_genes) %>%
  dplyr::mutate(Nanostring_gene_rank = row_number()) %>%
  dplyr::left_join(mouse_gene_ids, by = "genes") %>%
  dplyr::select(genes, ensembl_gene_id, external_synonym, everything()) %>%
  as.data.frame() -> DE_Nanostring_sig_genes

DE_Nanostring_sig_genes %>%
  dplyr::filter(is.na(ensembl_gene_id)) %>%
  dplyr::select(-c(external_synonym, ensembl_gene_id)) %>%
  dplyr::rename(external_synonym = genes) %>%
  dplyr::left_join(mouse_gene_ids, by = "external_synonym") %>%
  dplyr::select(genes, ensembl_gene_id, external_synonym, everything()) %>%
  as.data.frame() -> DE_Nanostring_sig_genes_NA

## read Gitler_GSE161621 results
readRDS("~/Documents/tmp/dacruz/202403_piol/data/new_DE/GSE161621_Blum_author_ann_DE_final.rds") %>% ## new DE
  dplyr::rename(genes = names(.)[1]) %>%
  dplyr::mutate(Gitler_GSE161621_gene_rank = row_number()) %>%
  dplyr::left_join(mouse_gene_ids, by = "genes") %>%
  dplyr::select(genes, ensembl_gene_id, external_synonym, everything()) %>%
  as.data.frame() -> Gitler_GSE161621_DE_RES

Gitler_GSE161621_DE_RES %>%
  dplyr::filter(is.na(ensembl_gene_id)) %>%
  dplyr::select(-c(external_synonym, ensembl_gene_id)) %>%
  dplyr::rename(external_synonym = genes) %>%
  dplyr::left_join(mouse_gene_ids, by = "external_synonym") %>%
  dplyr::select(genes, ensembl_gene_id, external_synonym, everything()) %>%
  as.data.frame() -> Gitler_GSE161621_DE_RES_NA

## read Le_Pichon_GSE167597 results
readRDS("~/Documents/tmp/dacruz/202403_piol/data/new_DE/GSE167597_Alkaslasi_author_ann_v4_DE.rds") %>% ## new CORRECT DE
  dplyr::rename(genes = names(.)[1]) %>%
  dplyr::mutate(Le_Pichon_GSE167597_gene_rank = row_number()) %>%
  dplyr::left_join(mouse_gene_ids, by = "genes") %>%
  dplyr::select(genes, ensembl_gene_id, external_synonym, everything()) %>%
  as.data.frame() -> Le_Pichon_GSE167597_DE_RES

Le_Pichon_GSE167597_DE_RES %>%
  dplyr::filter(is.na(ensembl_gene_id)) %>%
  dplyr::select(-c(external_synonym, ensembl_gene_id)) %>%
  dplyr::rename(external_synonym = genes) %>%
  dplyr::left_join(mouse_gene_ids, by = "external_synonym") %>%
  dplyr::select(genes, ensembl_gene_id, external_synonym, everything()) %>%
  as.data.frame() -> Le_Pichon_GSE167597_DE_RES_NA

## read Levine_GSE103892 results
readRDS(file = paste0(fp, "DE_results/RDS/GSE103892_DE_RES.rds")) %>%
  dplyr::rename(genes = names(.)[1]) %>%
  dplyr::mutate(Levine_GSE103892_gene_rank = row_number()) %>%
  dplyr::filter(group == "Chat_pos") %>%
  dplyr::left_join(mouse_gene_ids, by = "genes") %>%
  dplyr::select(genes, ensembl_gene_id, external_synonym, everything()) %>%
  as.data.frame() -> Levine_GSE103892_DE_RES

Levine_GSE103892_DE_RES %>%
  dplyr::filter(is.na(ensembl_gene_id)) %>%
  dplyr::select(-c(external_synonym, ensembl_gene_id)) %>%
  dplyr::rename(external_synonym = genes) %>%
  dplyr::left_join(mouse_gene_ids, by = "external_synonym") %>%
  dplyr::select(genes, ensembl_gene_id, external_synonym, everything()) %>%
  as.data.frame() -> Levine_GSE103892_DE_RES_NA

## de duplicate genes
dplyr::bind_rows(DE_Nanostring_sig_genes, DE_Nanostring_sig_genes_NA) %>%
  dplyr::filter(!is.na(ensembl_gene_id),
                !duplicated(ensembl_gene_id)) %>%
  as.data.frame() -> DE_Nanostring_ensembl

dplyr::bind_rows(Gitler_GSE161621_DE_RES, Gitler_GSE161621_DE_RES_NA) %>%
  dplyr::filter(!is.na(ensembl_gene_id),
                !duplicated(ensembl_gene_id)) %>%
  as.data.frame() -> Gitler_GSE161621_ensembl

dplyr::bind_rows(Le_Pichon_GSE167597_DE_RES, Le_Pichon_GSE167597_DE_RES_NA) %>%
  dplyr::filter(!is.na(ensembl_gene_id),
                !duplicated(ensembl_gene_id)) %>%
  as.data.frame() -> Le_Pichon_GSE167597_ensembl

dplyr::bind_rows(Levine_GSE103892_DE_RES, Levine_GSE103892_DE_RES_NA) %>%
  dplyr::filter(!is.na(ensembl_gene_id),
                !duplicated(ensembl_gene_id)) %>%
  as.data.frame() -> Levine_GSE103892_ensembl

## Venn diagram official gene symbols
list(`Piol`      = DE_Nanostring_ensembl$ensembl_gene_id,
     `Blum`      = Gitler_GSE161621_ensembl$ensembl_gene_id,
     `Alkaslasi` = Le_Pichon_GSE167597_ensembl$ensembl_gene_id,
     `Yadav`     = Levine_GSE190442_DE_RES2$ensembl_gene_id,
     `Gautier`   = Gautier_GSE228778_DE_RES2$ensembl_gene_id
     # Levine_GSE103892    = Levine_GSE103892_ensembl$ensembl_gene_id,
     ) %>%
  display_venn(
    main = paste0("Venn diagram of Chat+/CHAT+ vs Chat-/CHAT- DE genes mapped",
                  "to official mouse Ensembl gene IDs from Biomart"))

############################ DE genes unique to  ######################################
DE_Nanostring_ensembl %>%
    dplyr::filter(!ensembl_gene_id %in% Le_Pichon_GSE167597_ensembl$ensembl_gene_id, ## alkaslasi
                  !ensembl_gene_id %in% Gitler_GSE161621_ensembl$ensembl_gene_id, ## blum
                  !ensembl_gene_id %in% Gautier_GSE228778_DE_RES2$ensembl_gene_id, ## gautier
                  !ensembl_gene_id %in% Levine_GSE190442_DE_RES2$ensembl_gene_id ## yadav
                  ) %>%
    as.data.frame() -> DE_Nanostring_sig_genes_unique
dim(DE_Nanostring_sig_genes_unique)

############################ create Excel table ######################################

### create list of all DEG lists
list(DE_Nanostring_sig_genes_unique %>% dplyr::select(-contains("rank")), ## piol unique
     readRDS(file = paste0(fp, "DE_results/RDS/DE_Nanostring_sig_genes.rds")) %>%
       dplyr::select(-contains("rank")), ## piol
     Gautier_GSE228778_DE_RES %>% dplyr::select(-contains("rank")), ## gautier
     Levine_GSE190442_DE_RES %>% dplyr::select(-contains("rank")), ## yadav
     readRDS("~/Documents/tmp/dacruz/202403_piol/data/new_DE/GSE161621_Blum_author_ann_DE_final.rds"
             ) %>% ## new blum DE
       dplyr::rename(genes = names(.)[1]),
     readRDS("~/Documents/tmp/dacruz/202403_piol/data/new_DE/GSE167597_Alkaslasi_author_ann_v4_DE.rds"
             ) %>% ## new Alkaslasi DE
       dplyr::rename(genes = names(.)[1]),
     readRDS(file = paste0(fp, "DE_results/RDS/GSE103892_DE_RES.rds")) %>%
       dplyr::rename(genes = names(.)[1]) %>%
       dplyr::filter(group == "Chat_pos")) -> l0

### names of Excel sheets
names(l0) <- c("Novel Piol Chat+ vs Chat-",
               "Piol Chat+ vs Chat-",
               "Gautier CHAT+ vs CHAT-",
               "Yadav CHAT+ vs CHAT-",
               "Blum Chat+ vs Chat-",
               "Alkaslasi Chat+ vs Chat-",
               "Sathymurthy Chat+ vs Chat-")

## save Excel
saveRDS(object = l0, file = paste0(rp, "all_Chatpos_vs_Chatneg_unfiltered_DEG_lists_final_correct.rds"))
write_xlsx(x = l0, path = paste0(rp, "all_Chatpos_vs_Chatneg_unfiltered_DEG_lists_final_correct.xlsx"))

############################# correlation #####################################

## make subsets to join
Gitler_GSE161621_ensembl %>%
  dplyr::rename(Blum_lfc = avg_log2FC) %>% 
  dplyr::filter(!duplicated(genes)) %>% 
  dplyr::select(genes, Blum_lfc) -> a1

Le_Pichon_GSE167597_ensembl %>%
  dplyr::rename(Alkaslasi_lfc = avg_log2FC) %>%
  dplyr::filter(!duplicated(genes)) %>% 
  dplyr::select(genes, Alkaslasi_lfc) -> a2

Levine_GSE190442_DE_RES2 %>%
  dplyr::rename(Yadav_lfc = avg_log2FC) %>%
  dplyr::filter(!duplicated(genes)) %>% 
  dplyr::select(genes, Yadav_lfc) -> a3

Gautier_GSE228778_DE_RES2 %>%
  dplyr::rename(Gautier_lfc = avg_log2FC) %>%
  dplyr::filter(!duplicated(genes)) %>% 
  dplyr::select(genes, Gautier_lfc) -> a4

## join DE gene results
DE_Nanostring_ensembl %>%
  dplyr::filter(!duplicated(genes)) %>%
  dplyr::rename(Piol_lfc = log2FoldChange) %>%
  dplyr::arrange(desc(Piol_lfc)) %>%
  dplyr::left_join(a1, by = "genes") %>%
  dplyr::left_join(a2, by = "genes") %>%
  dplyr::left_join(a3, by = "genes") %>%
  dplyr::left_join(a4, by = "genes") %>%
  dplyr::select(genes, ensembl_gene_id, contains("_lfc")) %>%
  as.data.frame() -> de_ns

## write merged table to Excel
write_xlsx(x = de_ns,
           path = paste0(rp, "all_Chatpos_vs_Chatneg_deduplicated_LFC_final_correct.xlsx"))

### Correlation scatterplots
de_ns %>%
  ggplot(aes(x = Piol_lfc, y = Blum_lfc, label = genes)) +
  geom_point(shape = 21, fill = "mediumpurple1", color = "slateblue4", size = 1, stroke = 1) +
  xlab("Piol, et al. Log2 FC") + ylab("Blum, et al. 2021 Log2 FC") +
  geom_text_repel(force = 2) +
  stat_smooth(method = "lm", se = TRUE, formula = y ~ poly(x, 1, raw = TRUE),
              color = "black", fill = "gray") +
  ggpubr::stat_cor(method = "spearman") +
  theme_classic() +
  ggtitle("Piol, et al. vs Blum, et al. 2021 Log2 FC correlation") -> p1

de_ns %>%
  ggplot(aes(x = Piol_lfc, y = Alkaslasi_lfc, label = genes)) +
  geom_point(shape = 21, fill = "mediumpurple1", color = "slateblue4", size = 1, stroke = 1) +
  xlab("Piol, et al. Log2 FC") + ylab("Alkaslasi, et al. 2021 Log2 FC") +
  geom_text_repel() +
  stat_smooth(method = "lm", se = TRUE, formula = y ~ poly(x, 1, raw = TRUE),
              color = "black", fill = "gray") +
  ggpubr::stat_cor(method = "spearman") +
  theme_classic() +
  ggtitle("Piol, et al. vs Alkaslasi, et al. 2021 Log2 FC correlation") -> p2

de_ns %>%
  ggplot(aes(x = Piol_lfc, y = Yadav_lfc, label = genes)) +
  geom_point(shape = 21, fill = "mediumpurple1", color = "slateblue4", size = 1, stroke = 1) +
  xlab("Piol, et al. Log2 FC") + ylab("Yadav, et al. 2023") +
  geom_text_repel() +
  stat_smooth(method = "lm", se = TRUE, formula = y ~ poly(x, 1, raw = TRUE),
              color = "black", fill = "gray") +
  ggpubr::stat_cor(method = "spearman") +
  theme_classic() +
  ggtitle("Piol, et al. vs Yadav, et al. 2023 Log2 FC correlation") -> p3

de_ns %>%
  ggplot(aes(x = Piol_lfc, y = Gautier_lfc, label = genes)) +
  geom_point(shape = 21, fill = "mediumpurple1", color = "slateblue4", size = 1, stroke = 1) +
  xlab("Piol, et al. Log2 FC") + ylab("Gautier, et al. 2023") +
  geom_text_repel() +
  stat_smooth(method = "lm", se = TRUE, formula = y ~ poly(x, 1, raw = TRUE),
              color = "black", fill = "gray") +
  ggpubr::stat_cor(method = "spearman") +
  theme_classic() +
  ggtitle("Piol, et al. vs Gautier, et al. 2023 Log2 FC correlation") -> p4

## print PNGs of correlations
ggsave(filename = "~/Documents/tmp/dacruz/202403_piol/correlation_Piol_v_Blum.png",
       device = "png", plot = p1)
ggsave(filename = "~/Documents/tmp/dacruz/202403_piol/correlation_Piol_v_Alkaslasi.png",
       device = "png", plot = p2)
ggsave(filename = "~/Documents/tmp/dacruz/202403_piol/correlation_Piol_v_Yadav.png",
       device = "png", plot = p3)
ggsave(filename = "~/Documents/tmp/dacruz/202403_piol/correlation_Piol_v_Gautier.png",
       device = "png", plot = p4)