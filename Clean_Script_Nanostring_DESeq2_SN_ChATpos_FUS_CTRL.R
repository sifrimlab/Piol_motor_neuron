library("writexl")
library("DESeq2")
library("apeglm")
library("dplyr")
library("tidyr")
library("reshape2")
library("readxl")

# setwd("C:/Users/u0140675/OneDrive - KU Leuven/DESeq2_Analysis_Nanostring/Count_Matrices/Raw/SN/")
count_matrix_SN <- as.data.frame(read_excel("Count_Matrix_Raw_SN.xlsx"))
rownames(count_matrix_SN) <- count_matrix_SN$Genes


# setwd("C:/Users/u0140675/OneDrive - KU Leuven/DESeq2_Analysis_Nanostring/Count_Matrices/Metadata/SN")
metadata_SN <- as.data.frame(read_excel("Metadata_SN.xlsx"))

metadata_SN_ChATpos <- subset(metadata_SN, segment == "ChATpos")
rownames(metadata_SN_ChATpos) <- metadata_SN_ChATpos$sample

common_SN_ChATpos <- intersect(rownames(metadata_SN_ChATpos), colnames(count_matrix_SN))
count_matrix_SN_ChATpos <- count_matrix_SN[, common_SN_ChATpos]


# DESeq2.

metadata_SN_ChATpos <- metadata_SN_ChATpos[order(match(metadata_SN_ChATpos$sample, colnames(count_matrix_SN_ChATpos))), ]

dds <- DESeqDataSetFromMatrix(countData = count_matrix_SN_ChATpos,
                              colData = metadata_SN_ChATpos,
                              design= ~ class)

dds

dds <- DESeq(dds)


resultsNames(dds)

res <- results(dds, contrast = c("class", "FUS", "CTRL"))


df <- as.data.frame(res)
df$Genes <- row.names(df)

# setwd("C:/Users/u0140675/OneDrive - KU Leuven/")
# write_xlsx(df, "Results_SN_ChATpos_FUS_Vs_CTRL.xlsx")



