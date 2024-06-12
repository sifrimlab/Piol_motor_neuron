library("NanoStringNCTools")
library("GeomxTools")
library("GeoMxWorkflows")
library("ggplot2")
library("tidyverse")
library("DESeq2")
library("dplyr")
library("ggforce")
library("knitr")
library("tidyr")
library("writexl")

datadir <- file.path("../../..")
DCCFiles <- dir(file.path(datadir, "dccs"), pattern = ".dcc$",
                full.names = TRUE, recursive = TRUE)
PKCFiles <- dir(file.path(datadir, "pkcs"), pattern = ".pkc$",
                full.names = TRUE, recursive = TRUE)
SampleAnnotationFile <- dir(file.path(datadir, "annotation"), pattern = ".xlsx$",
                            full.names = TRUE, recursive = TRUE)
demoData <-
   readNanoStringGeoMxSet(dccFiles = DCCFiles,
                         pkcFiles = PKCFiles,
                         phenoDataFile = SampleAnnotationFile,
                         phenoDataSheet = "Template",
                         phenoDataDccColName = "Sample_ID",
                         protocolDataColNames = c("aoi", "roi"),
                         experimentDataColNames = c("panel"))

demoData


pkcs <- annotation(demoData)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))



#Exploratory Analysis.
demoData
assayData(demoData)[["exprs"]]


fData_Genes <- fData(demoData)
counts <- assayData(demoData)[["exprs"]]
colnames(counts) <- gsub("-", "_", colnames(counts))
class(counts)
rownames(counts)
colnames(counts)

counts2 <- as.data.frame(counts)
class(counts2)
rownames(counts2)
colnames(counts2)

metadata <- pData(demoData)
row.names(metadata)
names(metadata)
row.names(metadata) <- gsub("-", "_", row.names(metadata))
metadata <- data.frame(apply(metadata, c(1,2), str_replace_all, "-", "_"))
names(metadata) <- str_replace_all(names(metadata), "[.]", "_")
metadata$"sample" <- rownames(metadata)




# Split the Dataframe.

counts2_SN <- counts2[,1:24]
counts2_SC <- counts2[,25:48]
metadata_SN <- metadata[1:24,]
metadata_SC <- metadata[25:48,]

metadata_SN_ChATpos <- subset(metadata_SN, segment == "ChATpos")
metadata_SC_ChATpos <- subset(metadata_SC, segment == "ChATpos")

metadata_SN_ChATneg <- subset(metadata_SN, segment == "ChATneg")
metadata_SC_ChATneg <- subset(metadata_SC, segment == "ChATneg")


metadata_SC_Controls <- subset(metadata_SC, class == "CTRL")
metadata_SN_Controls <- subset(metadata_SN, class == "CTRL")
metadata_Controls <- subset(metadata, class == "CTRL")


metadata_SN_Controls_ChATpos <- subset(metadata_SN_Controls, segment == "ChATpos")
metadata_SN_Controls_ChATneg <- subset(metadata_SN_Controls, segment == "ChATneg")

metadata_ChATpos <- subset(metadata, segment == "ChATpos")
metadata_ChATpos_CTRL <- subset(metadata, segment == "ChATpos" & class == "CTRL")


# Add Genes.

counts2$RTS_ID <- rownames(counts2)
counts2 <- merge(x = counts2, y = fData_Genes[, c("RTS_ID", "TargetName")], by = "RTS_ID", all = TRUE)
counts2$Genes <- counts2$TargetName
counts2$TargetName <- NULL
counts2$RTS_ID <- NULL
counts2 <- subset(counts2, Genes != "Gm10406" & Genes != "LOC118568634" & Genes != "D830030K20Rik" & Genes != "NegProbe-WTX")
rownames(counts2) <- counts2$Genes
counts2$Genes <- NULL


counts2_SN$RTS_ID <- rownames(counts2_SN)
counts2_SN <- merge(x = counts2_SN, y = fData_Genes[, c("RTS_ID", "TargetName")], by = "RTS_ID", all = TRUE)
counts2_SN$Genes <- counts2_SN$TargetName
counts2_SN$TargetName <- NULL
counts2_SN$RTS_ID <- NULL
counts2_SN <- subset(counts2_SN, Genes != "Gm10406" & Genes != "LOC118568634" & Genes != "D830030K20Rik" & Genes != "NegProbe-WTX")
rownames(counts2_SN) <- counts2_SN$Genes
counts2_SN$Genes <- NULL


counts2_SC$RTS_ID <- rownames(counts2_SC)
counts2_SC <- merge(x = counts2_SC, y = fData_Genes[, c("RTS_ID", "TargetName")], by = "RTS_ID", all = TRUE)
counts2_SC$Genes <- counts2_SC$TargetName
counts2_SC$TargetName <- NULL
counts2_SC$RTS_ID <- NULL
counts2_SC <- subset(counts2_SC, Genes != "Gm10406" & Genes != "LOC118568634" & Genes != "D830030K20Rik" & Genes != "NegProbe-WTX")
rownames(counts2_SC) <- counts2_SC$Genes
counts2_SC$Genes <- NULL


# Split Dataframes for Controls.

common_samples_SC_Controls <- intersect(rownames(metadata_SC_Controls), colnames(counts2_SC))
counts2_SC_Controls <- counts2_SC[, common_samples_SC_Controls]

common_samples_SN_Controls <- intersect(rownames(metadata_SN_Controls), colnames(counts2_SN))
counts2_SN_Controls <- counts2_SN[, common_samples_SN_Controls]



common_samples_SN_Controls_ChATpos <- intersect(rownames(metadata_SN_Controls_ChATpos), colnames(counts2_SN_Controls))
counts2_SN_Controls_ChATpos <- counts2_SN_Controls[, common_samples_SN_Controls_ChATpos]

common_samples_SN_Controls_ChATneg <- intersect(rownames(metadata_SN_Controls_ChATneg), colnames(counts2_SN_Controls))
counts2_SN_Controls_ChATneg <- counts2_SN_Controls[, common_samples_SN_Controls_ChATneg]


# Split Dataframes for ChATpos or ChATneg.

common_samples_SC_ChATpos <- intersect(rownames(metadata_SC_ChATpos), colnames(counts2_SC))
counts2_SC_ChATpos <- counts2_SC[, common_samples_SC_ChATpos]

common_samples_SN_ChATpos <- intersect(rownames(metadata_SN_ChATpos), colnames(counts2_SN))
counts2_SN_ChATpos <- counts2_SN[, common_samples_SN_ChATpos]


common_samples_SC_ChATneg <- intersect(rownames(metadata_SC_ChATneg), colnames(counts2_SC))
counts2_SC_ChATneg <- counts2_SC[, common_samples_SC_ChATneg]

common_samples_SN_ChATneg <- intersect(rownames(metadata_SN_ChATneg), colnames(counts2_SN))
counts2_SN_ChATneg <- counts2_SN[, common_samples_SN_ChATneg]


common_samples_ChATpos <- intersect(rownames(metadata_ChATpos), colnames(counts2))
counts2_ChATpos <- counts2[, common_samples_ChATpos]

common_samples_ChATpos_CTRL <- intersect(rownames(metadata_ChATpos_CTRL), colnames(counts2))
counts2_ChATpos_CTRL <- counts2[, common_samples_ChATpos_CTRL]



# Split Dataframes for ChATpos or ChATneg and FUS-Controls.

common_samples_SC_ChATpos_FUS_Controls <- intersect(rownames(metadata_SC_ChATpos_FUS_Controls), colnames(counts2_SC_ChATpos))
counts2_SC_ChATpos_FUS_Controls <- counts2_SC_ChATpos[, common_samples_SC_ChATpos_FUS_Controls]

common_samples_SN_ChATpos_FUS_Controls <- intersect(rownames(metadata_SN_ChATpos_FUS_Controls), colnames(counts2_SN_ChATpos))
counts2_SN_ChATpos_FUS_Controls <- counts2_SN_ChATpos[, common_samples_SN_ChATpos_FUS_Controls]







#.......................................................................................

# Create Pivot Long Objects.
counts_long <- counts2 %>% 
  pivot_longer(!Genes, 
               names_to = "sample", 
               values_to = "cts")
counts_long <- full_join(counts_long, metadata, by = "sample")


counts_long_SN <- counts2_SN %>% 
  pivot_longer(!Genes, 
               names_to = "sample", 
               values_to = "cts")
counts_long_SN <- full_join(counts_long_SN, metadata_SN, by = "sample")


counts_long_SC <- counts2_SC %>% 
  pivot_longer(!Genes, 
               names_to = "sample", 
               values_to = "cts")
counts_long_SC <- full_join(counts_long_SC, metadata_SC, by = "sample")



# Save the Count Matrices.



setwd("../../../")

write_xlsx(counts_long,"Count_Matrix_Long_Whole_Dataset.xlsx")
write_xlsx(counts_long_SN,"Count_Matrix_Long_Schiatic_Nerve.xlsx")
write_xlsx(counts_long_SC,"Count_Matrix_Long_Spinal_Cord.xlsx")





#Plot the Distributions.

plot <- counts_long %>%
          ggplot(aes(cts, colour = sample)) + 
          geom_freqpoly(binwidth = 1) +
          xlim(0,40)
plot
path <- setwd("../../../")
ggsave(filename = "plot.png", plot = plot, device = "png", path = path, width = 12, height = 6)

#Plot the Distributions SN.

plot <- counts_long_SN %>%
  ggplot(aes(cts, colour = sample)) + 
  geom_freqpoly(binwidth = 1) +
  xlim(0,40)
plot
path <- setwd("../../../")
ggsave(filename = "plot.png", plot = plot, device = "png", path = path, width = 12, height = 6)

#Plot the Distributions SC.

plot <- counts_long_SC %>%
  ggplot(aes(cts, colour = sample)) + 
  geom_freqpoly(binwidth = 1) +
  xlim(0,40)
plot
path <- setwd("../../../")
ggsave(filename = "plot.png", plot = plot, device = "png", path = path, width = 12, height = 6)


# Plot Distributions by Slide (SN Samples).
plot <- counts_long_SN %>%
  ggplot(aes(cts, colour = `slide name`)) + 
  geom_freqpoly(binwidth = 1) +
  xlim(0,40)
plot
path <- setwd("../../../")
ggsave(filename = "plot.png", plot = plot, device = "png", path = path, width = 12, height = 6)


# Plot Distributions by Slide (SC Samples).
plot <- counts_long_SC %>%
ggplot(aes(cts, colour = `slide name`)) + 
  geom_freqpoly(binwidth = 1) +
  xlim(0,40)
plot
path <- setwd("../../../")
ggsave(filename = "plot.png", plot = plot, device = "png", path = path, width = 12, height = 6)

# Plot Distributions by Replicate_ID (SN Samples).
plot <- counts_long_SN %>%
  ggplot(aes(cts, colour = `replicate_id`)) + 
  geom_freqpoly(binwidth = 1) +
  xlim(0,40)
plot
path <- setwd("../../../")
ggsave(filename = "plot.png", plot = plot, device = "png", path = path, width = 12, height = 6)
