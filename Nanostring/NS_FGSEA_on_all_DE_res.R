## fgsea - this takes a list of DE results and finds enriched DE terms from the
## NS_SC_Chatpos_v_neg_20240604.Rmd workbook and generates lists of 
## significant DE GO and KEGG terms

## load libaries
library("dplyr")
library("GO.db")
library("hypeR")
library("tibble")
library("ggplot2")
library("stringr")
library("org.Mm.eg.db")
library("readxl")
library("writexl")

## workflow variables
proj_name <- "Spinal_cord_Chatpos_vs_Chatneg_control"
sp <- "Mus musculus"

## file path
# fp <- "~/Documents/tmp/dacruz/202406_piol/data/"
fp <- "/path/to/Piol/github/repo/data/"

## NS_SC_Chatpos_v_neg_20240604_final DE results
de_res <- "NanoString_spinal_cord_CTRL_DE_res_2024-06-04.xlsx"

## Load Pathway Gene Sets. Gene sets that are mapped to gene ontology (GO) and
## KEGG pathway terms are loaded from the `msigdb` database, for more info use:
## hypeR::msigdb_info()
GO_BP <- hypeR::msigdb_gsets(species = sp, category = "C5", subcategory = "BP")
GO_CC <- hypeR::msigdb_gsets(species = sp, category = "C5", subcategory = "CC")
GO_MF <- hypeR::msigdb_gsets(species = sp, category = "C5", subcategory = "MF")
KEGG <- hypeR::msigdb_gsets(species = sp, category = "C2", subcategory = "CP:KEGG")

## read DE results
readxl::read_excel(path = paste0(fp, de_res)) %>%
  dplyr::rename(gene_name = genes) %>%
  dplyr::filter(!is.na(entrez_id),
                !duplicated(entrez_id),
                !is.na(gene_name),
                !duplicated(gene_name)) %>%
  dplyr::arrange(pvalue) %>%
  dplyr::select(gene_name, pvalue) %>%
  as.data.frame() -> de_list

## set NA p-values to 0.99999999, so genes are not filtered out of ranked list
de_list$pvalue[is.na(de_list$pvalue)] <- 0.99999999

## convert DE genes to ranked list
de_list %>% tibble::deframe() -> rank_list

## GSEA GO BP 
fgsea::fgseaMultilevel(pathways = GO_BP$genesets,
                       stats = rank_list,
                       minSize = 1,
                       maxSize = Inf,
                       eps = 0,
                       nPermSimple = 1000) -> GO_BP_GSEA

## GSEA GO CC
fgsea::fgseaMultilevel(pathways = GO_CC$genesets,
                       stats = rank_list,
                       minSize = 1,
                       maxSize = Inf,
                       eps = 0,
                       nPermSimple = 1000) -> GO_CC_GSEA

## GSEA GO MF
fgsea::fgseaMultilevel(pathways = GO_MF$genesets,
                       stats = rank_list,
                       minSize = 1,
                       maxSize = Inf,
                       eps = 0,
                       nPermSimple = 1000) -> GO_MF_GSEA

## GSEA KEGG
fgsea::fgseaMultilevel(pathways = KEGG$genesets, 
                       stats = rank_list,
                       minSize = 1,
                       maxSize = Inf,
                       eps = 0,
                       nPermSimple = 1000) -> KEGG_GSEA

## arrange results of p-adjusted values and save results
GO_BP_GSEA %>%
  as.data.frame() %>%
  dplyr::arrange(padj) %>%
  dplyr::mutate(leadingEdge = gsub(",", ";", paste0(leadingEdge))) -> GO_BP_GSEA_res
GO_CC_GSEA %>%
  as.data.frame() %>%
  dplyr::arrange(padj) %>%
  dplyr::mutate(leadingEdge = gsub(",", ";", paste0(leadingEdge))) -> GO_CC_GSEA_res
GO_MF_GSEA %>%
  as.data.frame() %>%
  dplyr::arrange(padj) %>%
  dplyr::mutate(leadingEdge = gsub(",", ";", paste0(leadingEdge))) -> GO_MF_GSEA_res
KEGG_GSEA %>%
  as.data.frame() %>%
  dplyr::arrange(padj) %>%
  dplyr::mutate(leadingEdge = gsub(",", ";", paste0(leadingEdge))) -> KEGG_GSEA_res

## write results to Excel
f1 <- list(GO_BP_GSEA_res, GO_CC_GSEA_res, GO_MF_GSEA_res, KEGG_GSEA_res)
names(f1) <- c("GO_BP", "GO_CC", "GO_MF", "KEGG")
writexl::write_xlsx(x = f1, path = paste0(fp, proj_name, "_sig_GO_KEGG_terms",
                                          "_", Sys.Date(), ".xlsx"))

################################################################################

## workflow variables
# proj_name <- "SN_NS_over_5_counts"
proj_name <- "Sciatic_nerve_Chatpos_vs_Chatneg_control"
sp <- "Mus musculus"

## file path
# fp <- "~/Documents/tmp/dacruz/202406_piol/data/"
fp <- "/path/to/Piol/github/repo/data/"

## NS_SC_Chatpos_v_neg_20240604_final DE results
de_res <- "NanoString_sciatic_nerve_CTRL_segment_DE_res_2024-06-04.xlsx"

## Load Pathway Gene Sets. Gene sets that are mapped to gene ontology (GO) and
## KEGG pathway terms are loaded from the `msigdb` database, for more info use:
## hypeR::msigdb_info()
GO_BP <- hypeR::msigdb_gsets(species = sp, category = "C5", subcategory = "BP")
GO_CC <- hypeR::msigdb_gsets(species = sp, category = "C5", subcategory = "CC")
GO_MF <- hypeR::msigdb_gsets(species = sp, category = "C5", subcategory = "MF")
KEGG <- hypeR::msigdb_gsets(species = sp, category = "C2", subcategory = "CP:KEGG")

## read DE results
readxl::read_excel(path = paste0(fp, de_res)) %>%
  dplyr::rename(gene_name = genes) %>%
  dplyr::filter(!is.na(entrez_id),
                !duplicated(entrez_id),
                !is.na(gene_name),
                !duplicated(gene_name)) %>%
  dplyr::select(gene_name, contains("dcc")) %>%
  tibble::column_to_rownames(var = "gene_name") %>%
  dplyr::mutate(mean_norm_counts = rowMeans(.)) %>%
  dplyr::arrange(desc(mean_norm_counts)) %>%
  tibble::rownames_to_column(var = "gene_name") %>%
  dplyr::select(gene_name, mean_norm_counts) %>%
  dplyr::filter(mean_norm_counts > 5) %>% 
  tibble::deframe() -> rank_list

## NOTE: so unlike the spinal cord GO DE script, we will be ranking the genes
## mean row-wise expression of the normalize DESeq counts (i.e. mean_norm_counts)
## The reason that we are doing this, is because we want to profile the GO terms
## of genes enriched in axons (i.e. sciatic nerve) and we want to effectively
## normalize by area, and so the DESeq2 size factor normalization is an effective
## proxy for normalizing expression over area, since we also cannot assume that
## there is even expression over the sciatic nerve GeoMx samples, and therefore
## normalizing by total sample gene count (i.e. size factors) is more appropriate

## GSEA GO BP
fgsea::fgseaMultilevel(pathways = GO_BP$genesets,
                       stats = rank_list,
                       minSize = 1,
                       maxSize = Inf,
                       eps = 0,
                       nPermSimple = 1000) -> GO_BP_GSEA

## GSEA GO CC
fgsea::fgseaMultilevel(pathways = GO_CC$genesets,
                       stats = rank_list,
                       minSize = 1,
                       maxSize = Inf,
                       eps = 0,
                       nPermSimple = 1000) -> GO_CC_GSEA

## GSEA GO MF
fgsea::fgseaMultilevel(pathways = GO_MF$genesets,
                       stats = rank_list,
                       minSize = 1,
                       maxSize = Inf,
                       eps = 0,
                       nPermSimple = 1000) -> GO_MF_GSEA

## GSEA KEGG
fgsea::fgseaMultilevel(pathways = KEGG$genesets,
                       stats = rank_list,
                       minSize = 1,
                       maxSize = Inf,
                       eps = 0,
                       nPermSimple = 1000) -> KEGG_GSEA

## save results
GO_BP_GSEA %>%
  as.data.frame() %>%
  dplyr::arrange(padj) %>%
  dplyr::mutate(
    leadingEdge = gsub(",", ";", paste0(leadingEdge))) -> GO_BP_GSEA_res
GO_CC_GSEA %>%
  as.data.frame() %>%
  dplyr::arrange(padj) %>%
  dplyr::mutate(
    leadingEdge = gsub(",", ";", paste0(leadingEdge))) -> GO_CC_GSEA_res
GO_MF_GSEA %>%
  as.data.frame() %>%
  dplyr::arrange(padj) %>%
  dplyr::mutate(
    leadingEdge = gsub(",", ";", paste0(leadingEdge))) -> GO_MF_GSEA_res
KEGG_GSEA %>%
  as.data.frame() %>%
  dplyr::arrange(padj) %>%
  dplyr::mutate(
    leadingEdge = gsub(",", ";", paste0(leadingEdge))) -> KEGG_GSEA_res

## write results to Excel
f1 <- list(GO_BP_GSEA_res, GO_CC_GSEA_res, GO_MF_GSEA_res, KEGG_GSEA_res)
names(f1) <- c("GO_BP", "GO_CC", "GO_MF", "KEGG")
writexl::write_xlsx(x = f1, path = paste0(fp, proj_name, "_sig_GO_KEGG_terms",
                                          "_", Sys.Date(), ".xlsx"))