library("DESeq2")
library("ggplot2")

# SC..............................................................................

setwd("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString_DaCruz_Alejandro/Maria/Analysis/Results/DESeq2_Analysis/Shrunken/SC/Excel_Files/No_Interaction_Terms/")
SC_Chat <- read_excel("Results_ChATpos_Vs_ChATneg_Shrunk_Ordered.xlsx")


#SC_Chat_Sub <- subset(SC_Chat, padj < 0.05 & (log2FoldChange > 3.2 | log2FoldChange < -3.2))


setwd("C:/Users/u0140675/OneDrive - KU Leuven/DESeq2_Analysis/Plots/Normalized_Count_ScatterPlots/SC/All_Genes/")

# Define a vector of gene IDs or names to plot
gene_list <- SC_Chat$Genes

# Define the two grouping variables
group1 <- "segment"
group2 <- "class"

# Create plots for each gene.


tic()

for (gene in gene_list) {
  # Subset the metadata to include only the two grouping variables and the current gene
  metadata_SC_sub <- metadata_SC[,c(group1, group2)]
  metadata_SC_sub$gene <- gene
  
  # Merge the normalized counts with the subsetted metadata
  df <- data.frame(counts(dds, normalized = TRUE)[gene,], metadata_SC_sub)
  colnames(df) <- c("norm_counts", group1, group2, "gene")
  
  # Create a scatterplot of the normalized counts using ggplot
  p <- ggplot(df, aes(x = reorder(factor(df[,group1]), df[,group1], median), y = norm_counts, shape = factor(df[,group1]), color = factor(df[,group2]))) +
    geom_point(size = 2) +
    facet_grid(~ factor(df[,group2]), scales = "free") +
    labs(title = gene, x = group1, y = "Normalized Counts", shape = group1, color = group2)
  
  # Save the plot as a PDF file
  ggsave(paste(gene, ".png", sep = ""), p)
}


toc()

# SN..............................................................................


setwd("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString_DaCruz_Alejandro/Maria/Analysis/Results/DESeq2_Analysis/Shrunken/SN/Excel_Files/No_Interaction_Terms/")
SN_Chat <- read_excel("Results_ChATpos_Vs_ChATneg_Shrunk_Ordered.xlsx")


#SN_Chat_Sub <- subset(SN_Chat, padj < 0.05 & (log2FoldChange > 3.2 | log2FoldChange < -3.2))


setwd("C:/Users/u0140675/OneDrive - KU Leuven/DESeq2_Analysis/Plots/Normalized_Count_ScatterPlots/SN/All_Genes/")

# Define a vector of gene IDs or names to plot
gene_list <- SN_Chat$Genes

# Define the two grouping variables
group1 <- "segment"
group2 <- "class"

# Create plots for each gene

tic()

for (gene in gene_list) {
  # Subset the metadata to include only the two grouping variables and the current gene
  metadata_SN_sub <- metadata_SN[,c(group1, group2)]
  metadata_SN_sub$gene <- gene
  
  # Merge the normalized counts with the subsetted metadata
  df <- data.frame(counts(dds, normalized = TRUE)[gene,], metadata_SN_sub)
  colnames(df) <- c("norm_counts", group1, group2, "gene")
  
  # Create a scatterplot of the normalized counts using ggplot
  p <- ggplot(df, aes(x = reorder(factor(df[,group1]), df[,group1], median), y = norm_counts, shape = factor(df[,group1]), color = factor(df[,group2]))) +
    geom_point(size = 2) +
    facet_grid(~ factor(df[,group2]), scales = "free") +
    labs(title = gene, x = group1, y = "Normalized Counts", shape = group1, color = group2)
  
  # Save the plot as a PDF file
  ggsave(paste(gene, ".png", sep = ""), p)
}

toc()


# SC Vs SN........................................................................

setwd("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString_DaCruz_Alejandro/Maria/Analysis/Results/DESeq2_Analysis/Shrunken//Whole_Dataset/Excel_Files/")
data <- read_excel("Results_SN_Vs_SC_Shrunk_Ordered.xlsx")


#SN_Chat_Sub <- subset(SN_Chat, padj < 0.05 & (log2FoldChange > 3.2 | log2FoldChange < -3.2))


setwd("C:/Users/u0140675/OneDrive - KU Leuven/DESeq2_Analysis/Plots/Normalized_Count_ScatterPlots/SN_and_SC_(Used_for_the_DE_SC_Vs_SN)/All_Genes/")

# Define a vector of gene IDs or names to plot
gene_list <- data$Genes

# Define the two grouping variables
group1 <- "segment"
group2 <- "class"
group3 <- "region"

# Create plots for each gene

for (gene in gene_list) {
  # Subset the metadata to include all three grouping variables and the current gene
  metadata_sub <- metadata[,c("segment", "class", "region")]
  metadata_sub$gene <- gene
  
  # Merge the normalized counts with the subsetted metadata
  df <- data.frame(counts(dds, normalized = TRUE)[gene,], metadata_sub)
  colnames(df) <- c("norm_counts", "segment", "class", "region", "gene")
  
  # Create a scatterplot of the normalized counts using ggplot
  p <- ggplot(df, aes(x = interaction(class, region, lex.order = TRUE), y = norm_counts, shape = segment, color = class)) +
    geom_point(size = 2) +
    scale_y_log10() +
    facet_grid(~region*segment, scales = "free_x") +
    labs(title = gene, x = "Class per Region", y = "log10 Normalized Counts", shape = "Segment", color = "Class")
  
  # Save the plot as a PNG file
  ggsave(paste(gene, ".png", sep = ""), p, width = 10, height = 6)
}
