library("pheatmap")
library("gplots")
library("ggplot2")
library("dplyr")
library("ggplotify")
library("gridExtra")
library("readxl")
library("RColorBrewer")
library("stats")



setwd("../../../")
df <- read_excel("DEG_List.xlsx")
df <- subset(df, padj < 0.05)
df <- df[order(-df$log2FoldChange),]



df_top_chat_pos <- df[1:20,]
df_top_chat_neg <- tail(df, 20)
df_top <- merge(df_top_chat_pos, df_top_chat_neg, all = TRUE)
df_top <- df_top[order(-df_top$log2FoldChange),]


genes <- df_top$genes
genes



counts <- counts(dds, normalized = TRUE) ##Use dds from DESeq.

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

gene_order <- rev(unique(df_top$genes))
gene_order <- unique(df_top$genes)



heatmap <- ggplot(counts_df_sub_long, aes(x = factor(sample), y = factor(Genes, levels = gene_order), fill = cts)) +
  geom_tile(colour="white", size = 0.5) +
  facet_grid(~ segment, scales = "free") +
  scale_x_discrete(labels = NULL) +
  scale_fill_gradientn(colours = colorRampPalette(c("white", "slateblue4","khaki1", "violetred4"))(100)) +
  labs(x = "Samples",
       y = "Genes",
       title = "Heatmap of the Top 20 Differentialy Expressed Genes",
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




setwd("../../../")
ggsave("Heatmap.png", plot = heatmap, width = 10, height = 6)

