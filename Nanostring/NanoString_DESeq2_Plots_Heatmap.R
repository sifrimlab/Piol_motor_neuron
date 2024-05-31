library("pheatmap")
library("gplots")
library("ggplot2")
library("dplyr")
library("ggplotify")
library("gridExtra")
library("readxl")
library("RColorBrewer")
library("stats")


setwd("C:/Users/u0140675/OneDrive - KU Leuven/DESeq2_Analysis_Nanostring/No_Shrunken/SC/Excel_Files/With_Interaction_Terms/CTRL")
df <- read_excel("Results_SC_ChATpos_Vs_ChATneg_CTRL_Ordered.xlsx")
df <- subset(df, pvalue < 0.05)
df <- df[order(-df$log2FoldChange),]

df_top_chat_pos <- df[1:50,]
df_top_chat_neg <- tail(df, 50)
df_top <- merge(df_top_chat_pos, df_top_chat_neg, all = TRUE)
df_top <- df_top[order(-df_top$log2FoldChange),]


genes <- df_top_chat_neg$Genes
genes <- df_top$Genes
genes

vst <- vst(dds, blind = FALSE)
counts_df <- assay(vst) 

counts_df_sub <- counts_df[rownames(counts_df) %in% genes, ]
counts_df_sub <- counts_df_sub[df_top$Genes, ]
counts_df_sub_log <- as.data.frame(log2(counts_df_sub))
counts_df_sub_log$Genes <- rownames(counts_df_sub_log)


counts_df_sub_log_long <- counts_df_sub_log %>% 
  pivot_longer(!Genes, 
               names_to = "sample", 
               values_to = "cts")
counts_df_sub_log_long <- full_join(counts_df_sub_log_long, metadata_SN_Controls, by = "sample")

gene_order <- rev(unique(counts_df_sub_log_long$Genes))


counts_df_sub_log_long$segment_2 <- paste0(counts_df_sub_log_long$segment, " Segments")
counts_df_sub_log_long$segment_2 <- reorder(counts_df_sub_log_long$segment_2, 
                                          match(counts_df_sub_log_long$segment, 
                                                metadata_SN_Controls$segment))


# or.................................................................................................

counts <- counts(dds, normalized = TRUE)

# Perform a log transformation (e.g., log2) on the count matrix

log_count_matrix <- log2(counts + 1)

# Calculate the mean count value for each gene
mean_counts <- rowMeans(log_count_matrix)

# Center the counts around zero by subtracting the mean count value for each gene
counts_df <- sweep(log_count_matrix, 1, mean_counts, "-")

counts_df_sub <- as.data.frame(counts_df[rownames(counts_df) %in% genes, ])
counts_df_sub$Genes <- rownames(counts_df_sub)



counts_df_sub_long <- counts_df_sub %>% 
  pivot_longer(!Genes, 
               names_to = "sample", 
               values_to = "cts")
counts_df_sub_long <- full_join(counts_df_sub_long, metadata_SC_Controls, by = "sample")

gene_order <- rev(unique(genes))
gene_order <- unique(genes)



heatmap <- ggplot(counts_df_sub_long, aes(x = factor(sample), y = factor(Genes, levels = gene_order), fill = cts)) +
  geom_tile(colour="white", size = 0.5) +
  facet_grid(~ segment, scales = "free") +
  scale_x_discrete(labels = NULL) +
  scale_fill_gradientn(colours = colorRampPalette(c("white", "slateblue4","khaki1", "violetred4"))(100)) +
  #scale_fill_gradientn(colours = colorRampPalette(c("white", "violetred4"))(100)) +
  #scale_fill_gradientn(colours = colorRampPalette(c("white", "red"))(100)) +
  #scale_fill_gradientn(colours = colorRampPalette(c("slateblue4", "white", "violetred4"))(100)) +
  labs(x = "Samples",
       y = "Genes",
       title = "Heatmap of the Top 50 Differentialy Expressed Genes",
       subtitle = "Spinal Cord - Controls - ChATneg Vs ChATpos",
       fill = "log2 Normalized Counts") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.2, size = 12),
        axis.text.y = element_text(size = 9.5),
        axis.title.x = element_text(color = "gray18"),
        axis.title.y = element_text(color = "gray18"), 
        plot.title = element_text(color = "gray18", hjust = 0.1, size = 13),
        plot.subtitle = element_text(color = "gray18", hjust = 0.1, size = 10),
        axis.line.x = element_blank(),
        axis.ticks.y = element_line(color = "midnightblue", linewidth = 0.4),
        axis.ticks.x = element_blank())

heatmap


  
setwd("C:/Users/u0140675/OneDrive - KU Leuven/DESeq2_Analysis_Nanostring/Plots/Heatmaps/Without_Shrinkage/Color_1")
ggsave("Heatmap_SC_Controls_ChATpos_Vs_ChATneg_top_50_ChATneg.png", plot = heatmap, width = 10, height = 9)




# Heatmap of Specific Genes.

setwd("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString_DaCruz_Alejandro/Maria/My_Stuff/")
df <- read_excel("Genes_Selected_by_Diana_SC_ChATpos_Vs_ChATneg.xlsx")

genes <- df$Genes


vst <- vst(dds, blind = FALSE)
counts_df <- assay(vst) 


counts_df_sub <- counts_df[rownames(counts_df) %in% genes, ]
counts_df_sub <- counts_df_sub[df$Genes, ]
counts_df_sub_log <- as.data.frame(log2(counts_df_sub))
counts_df_sub_log$Genes <- rownames(counts_df_sub_log)



counts_df_sub_log_long <- counts_df_sub_log %>% 
  pivot_longer(!Genes, 
               names_to = "sample", 
               values_to = "cts")
counts_df_sub_log_long <- full_join(counts_df_sub_log_long, metadata_SC, by = "sample")

gene_order <- rev(unique(counts_df_sub_log_long$Genes))


counts_df_sub_log_long$segment_2 <- paste0(counts_df_sub_log_long$segment, " Segments")
counts_df_sub_log_long$segment_2 <- reorder(counts_df_sub_log_long$segment_2, 
                                            match(counts_df_sub_log_long$segment, 
                                                  metadata_SC$segment))


class_order <- c("TDP_43", "CTRL", "FUS")
counts_df_sub_log_long$class <- factor(counts_df_sub_log_long$class, levels = class_order)
sorted_df <- counts_df_sub_log_long[order(counts_df_sub_log_long$Genes, counts_df_sub_log_long$class), ]


sorted_df$sample <- factor(sorted_df$sample, levels = unique(sorted_df$sample))


heatmap <- ggplot(sorted_df, aes(x = sample, 
                                 y = factor(Genes, levels = gene_order), 
                                 fill = cts)) +
  geom_tile(colour="white", size = 0.5) +
  facet_grid(~ segment_2, scales = "free") +
  scale_fill_gradientn(colours = colorRampPalette(c("white", "slateblue4", "khaki1", "violetred4"))(100)) +
  #scale_fill_gradientn(colours = colorRampPalette(c("white", "violetred4"))(100)) +
  #scale_fill_gradientn(colours = colorRampPalette(c("white", "red"))(100)) +
  #scale_fill_gradientn(colours = colorRampPalette(c("slateblue4", "white", "violetred4"))(100)) +
  labs(x = "Samples",
       y = "Genes",
       title = "Heatmap of Differentialy Expressed Genes",
       subtitle = "Spinal Cord -  ChATpos Vs ChATneg",
       fill = "log2 Normalized Counts") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.2, size = 12),
        axis.text.y = element_text(size = 9.5),
        axis.title.x = element_text(color = "gray18"),
        axis.title.y = element_text(color = "gray18"), 
        plot.title = element_text(color = "gray18", hjust = 0.1, size = 13),
        plot.subtitle = element_text(color = "gray18", hjust = 0.1, size = 10),
        axis.line.x = element_blank(),
        axis.ticks.y = element_line(color = "midnightblue", linewidth = 0.4),
        axis.ticks.x = element_blank()) +
        #axis.ticks.x = element_line(color = "midnightblue", linewidth = 0.4)) +
  scale_x_discrete(breaks = function(x) x[seq(2, length(x), by = 4)],
                   labels = c("TDP-43", "Control", "FUS", "TDP-43", "Control", "FUS")) +
  geom_vline(xintercept = c(4.5, 8.5), color = "white", size = 3)
  
heatmap


setwd("C:/Users/u0140675/OneDrive - KU Leuven/DESeq2_Analysis_Nanostring/Plots/Heatmaps/Without_Shrinkage/Color_1")
ggsave("Heatmap_SC_ChATpos_Vs_ChATneg_Selected_Genes_ChATneg.png", plot = heatmap, width = 10, height = 6)


#..........................................................................................................





