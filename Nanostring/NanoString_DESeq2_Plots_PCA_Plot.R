library("readxl")
library("ggplot2")
library("ggplotify")
library("ggrepel")
library("ggforce")
library("ggalt")
library("dplyr")
library("cowplot")
library("stats")

# PCA Plots. SC - SN.........................................................................

# Perform variance stabilizing transformation
vst <- vst(dds, blind = FALSE)

# Create PCA plot.
pcaData <- plotPCA(vst, intgroup = c("region", "class", "segment"), returnData = TRUE)


pcaData_grouped <- pcaData %>%
  group_by(region)

# Create ggplot object with shape and color aesthetics.

setwd("C:/Users/u0140675/OneDrive - KU Leuven/DESeq2_Analysis_Nanostring/Plots/PCA_Plots/PCA_Plots_Controls_Fus/")

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

ggsave(pcaPlot, filename = "PCA1_PCA2_SN_SC_FUS_Controls.png", width = 8, height = 6)


#...........................................................................................


# PCA Plots. SC or SN........................................................................

# Perform variance stabilizing transformation
vst <- vst(dds, blind = FALSE)

# Create PCA plot.

plotPCA(vst, intgroup = c("class", "segment"))


# Create ggplot object with shape and color aesthetics.

setwd("C:/Users/u0140675/OneDrive - KU Leuven/DESeq2_Analysis_Nanostring/Plots/PCA_Plots/PCA_Plots_Controls_Fus/")

percentVar <- round(100 * attr(pcaData, "percentVar"))



pcaData <- plotPCA(vst, intgroup = c("class", "segment"), returnData = TRUE)

pcaPlot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = class, shape = segment)) + 
  geom_point(size = 3) + 
  theme_bw() +
  labs(x = paste0("PC1: ",percentVar[1],"% variance"), 
       y = paste0("PC2: ",percentVar[2],"% variance"),
       title = "PCA plot",
       subtitle = "Visualization of the First Two Principal Components (Sciatic Nerve)",
       color = "Class", shape = "Segment")

pcaPlot

ggsave(pcaPlot, filename = "PCA1_PCA2_SN_FUS_Controls.png", width = 8, height = 6)


# PCA Plots. SC or SN and ChATpos or ChATneg........................................................................

# Perform variance stabilizing transformation.

vst <- vst(dds, blind = FALSE)

# Create PCA plot.

plotPCA(vst, intgroup = c("class"))

pcaData <- plotPCA(vst, intgroup = c("class"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Create ggplot object with shape and color aesthetics.

setwd("C:/Users/u0140675/OneDrive - KU Leuven/DESeq2_Analysis_Nanostring/Plots/PCA_Plots/PCA_Plots_Controls_Fus/")


pcaPlot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = class)) + 
  geom_point(size = 3) + 
  theme_bw() +
  labs(x = paste0("PC1: ",percentVar[1],"% variance"), 
       y = paste0("PC2: ",percentVar[2],"% variance"),
       title = "PCA plot",
       subtitle = "Visualization of the First Two Principal Components (Sciatic Nerve ChATneg)",
       color = "Class")

pcaPlot

ggsave(pcaPlot, filename = "PCA1_PCA2_SN_ChATneg_FUS_Controls.png", width = 8, height = 6)


# Plot more PCAs.....................................................................................................
## Maybe Wrong.

vst <- vst(dds, blind = FALSE)
assay(vst)
vst_df <- assay(vst)
vst_df
pca <- prcomp(t(vst_df))
summary(pca)
pca$x

percentVar_2 <- round(100 * (cumsum(pca$sdev^2) / sum(pca$sdev^2)))

pca_df <- cbind(metadata_SN, pca$x)

pcaPlot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = class, shape = segment)) + 
  geom_point(size = 3) + 
  theme_bw() +
  labs(x = paste0("PC1: ",percentVar_2[1],"% variance"), 
       y = paste0("PC2: ",percentVar_2[2],"% variance"),
       title = "PCA plot",
       subtitle = "Visualization of the Principal Components 3 and 4(Sciatic Nerve)",
       color = "Class", shape = "Segment")

pcaPlot


# Create data frame with metadata and PC3 and PC4 values for input to ggplot
df <- cbind(meta, pca$x)
ggplot(df) + geom_point(aes(x=PC3, y=PC4, color = sampletype))

