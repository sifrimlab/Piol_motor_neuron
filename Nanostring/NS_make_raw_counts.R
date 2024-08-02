## this script takes the long for Nanostring GeoMx counts processed by Maria in
## her script(s): to be assigned
## ... and coverts this table into a "wide" count matrix where the rows are genes
## and the columns correspond to the Nanostring GeoMx 

## load libaries
library("dplyr")
library("tidyr")
library("tibble")
library("readxl")
library("readr")

## input file path
# fp <- "~/Documents/tmp/dacruz/202311_piol/data/"
fp <- "/path/to/Piol/github/repo/data/"

## output file path
# out_path <- "~/Documents/tmp/dacruz/202406_piol/data/"
out_path <- "/path/to/Piol/github/repo/data/"

## read long format count file
readxl::read_excel(path = paste0(fp, "count_matrix_Long_Whole_Dataset.xlsx")) %>% 
  as.data.frame() -> raw_df

## convert long data to metadata table (colData) for DESeq2
raw_df %>% 
  dplyr::filter(!duplicated(sample)) %>%
  dplyr::select(-c(Genes, cts, expected_neg, scan_name)) %>%
  dplyr::filter(class %in% c("CTRL", "FUS")) %>%
  dplyr::mutate_all(as.factor) %>%
  as.data.frame() -> meta_data
rownames(meta_data) <- meta_data$sample

## convert long format dataframe to wide count matrix
raw_df %>%
  dplyr::select(Genes, sample, cts) %>%
  tidyr::pivot_wider(names_from = "Genes", values_from = "cts") %>%
  tibble::column_to_rownames(var = "sample") %>%
  t() %>% as.data.frame() %>%
  dplyr::select(meta_data$sample) %>%
  as.data.frame() -> cts

## save metadata and processed counts as TSV files
write.table(x = meta_data, file = paste0(out_path, "Nanostring_metadata_table.tsv"),
            sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
write.table(x = cts, file = paste0(out_path, "Nanostring_processed_count_matrix.tsv"),
            sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)

## visual check of files
# test_meta <- read.table(file = paste0(out_path, "Nanostring_metadata_table.tsv"),
#                         header = TRUE, sep = "\t")
# View(test_meta)