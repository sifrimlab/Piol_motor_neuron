library("clusterProfiler")
library("pathview")
library("ggplot2")
library("dplyr")
library("readxl")
library("stats")

organism = "org.Mm.eg.db"
library(organism, character.only = TRUE)



setwd("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString_DaCruz_Alejandro/Maria/Analysis/Results/DESeq2_Analysis/Shrunken/SN/Excel_Files/With_Interaction_Terms/ChATpos")
df <- read_excel("Results_SN_TDP_43_Vs_FUS_ChATpos_Shrunk_Ordered.xlsx")

original_gene_list <- df$log2FoldChange
names(original_gene_list) <- df$Genes
gene_list <- na.omit(original_gene_list)
gene_list <- sort(gene_list, decreasing = TRUE)


keytypes(org.Mm.eg.db)



gse <- gseGO(geneList = gene_list, 
                        ont = "ALL", 
                        keyType = "SYMBOL",
                        minGSSize = 3, 
                        maxGSSize = 800, 
                        pvalueCutoff = 0.05, 
                        verbose = TRUE, 
                        OrgDb = organism, 
                        pAdjustMethod = "none")

warnings()


require(DOSE)

plot <- dotplot(gse, showCategory = 20, split = ".sign") + facet_grid(.~.sign)
path <- setwd("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString_DaCruz_Alejandro/Maria/Analysis/Results/Functional_Analysis/Shrunken/SN/DotPlots/With_Interaction_Terms/ChATpos")

ggsave(filename = "DotPlot_TDP_43_Vs_FUS_ChATpos_SN.png", plot = plot, device = "png", path = path, width = 10, height = 28)



