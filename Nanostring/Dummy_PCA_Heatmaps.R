# Load the RNA-seq data into R as a matrix or data frame
rna_data <- read.table("rna_data.txt", header = TRUE, row.names = 1)

# Normalize the data using ex. TMM normalization. (Or RLE, Quantile or whatever)
library(edgeR)
cpm_data <- cpm(rna_data)
norm_data <- calcNormFactors(cpm_data)

# Transpose the data. So samples in columns and genes in rows.
trans_data <- t(norm_data)

#...................................................................................................

# PCA.

# Perform PCA
pca_data <- prcomp(trans_data)

# Create a plot of the first two principal components
library(ggplot2)
pc_df <- data.frame(PC1=pca_data$x[,1], PC2=pca_data$x[,2], group=colnames(trans_data))
ggplot(pc_df, aes(x=PC1, y=PC2, color=group)) + 
  geom_point(size=3) + 
  labs(x="PC1", y="PC2", title="PCA Plot")


# Example of ggplot2 editing.

percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaPlot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = class, shape = segment, group = region)) + 
  geom_point(size = 3) + 
  theme_bw() +
  labs(x = paste0("PC1: ",percentVar[1],"% variance"), 
       y = paste0("PC2: ",percentVar[2],"% variance"),
       title = "PCA plot",
       subtitle = "Visualization of the First Two Principal Components",
       color = "Class", shape = "Segment") +
  annotate("text", x = 12, y = 7.8, label = "Spinal Cord", color = "blue") +
  annotate("text", x = 12, y = -7.8, label = "Spinal Cord", color = "blue") +
  annotate("text", x = -17.5, y = 3.5, label = "Sciatic Nerve", color = "purple")

pcaPlot

ggsave(pcaPlot, filename = "PCA1_PCA2_SN_SC.png", width = 8, height = 6)

#.....................................................................................................


# Heatmaps. 

pheatmap(trans_data, scale = "row", color = "bluewhiteorange", clustering_distance_rows = "correlation", clustering_distance_cols = "correlation")

# annotation_row = .... , annotation_column = ......