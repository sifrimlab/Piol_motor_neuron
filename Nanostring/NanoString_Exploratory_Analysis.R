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



datadir <- file.path("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString_DaCruz_Alejandro/Maria/Files/WTA")
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

str(demoData)

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

metadata_SC_FUS_Controls <- subset(metadata_SC, class != "TDP_43")
metadata_SN_FUS_Controls <- subset(metadata_SN, class != "TDP_43")
metadata_FUS_Controls <- subset(metadata, class != "TDP_43")

metadata_SC_Controls <- subset(metadata_SC, class == "CTRL")
metadata_SN_Controls <- subset(metadata_SN, class == "CTRL")
metadata_Controls <- subset(metadata, class == "CTRL")


metadata_SN_Controls_ChATpos <- subset(metadata_SN_Controls, segment == "ChATpos")
metadata_SN_Controls_ChATneg <- subset(metadata_SN_Controls, segment == "ChATneg")

metadata_ChATpos <- subset(metadata, segment == "ChATpos")
metadata_ChATpos_CTRL <- subset(metadata, segment == "ChATpos" & class == "CTRL")

metadata_SN_ChATpos_FUS_Controls <- subset(metadata_SN_ChATpos, class != "TDP_43")
metadata_SC_ChATpos_FUS_Controls <- subset(metadata_SC_ChATpos, class != "TDP_43")

metadata_SN_ChATneg_FUS_Controls <- subset(metadata_SN_ChATneg, class != "TDP_43")
metadata_SC_ChATneg_FUS_Controls <- subset(metadata_SC_ChATneg, class != "TDP_43")



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

# Split Dataframes for FUS and Controls.

common_samples_SC_FUS_Controls <- intersect(rownames(metadata_SC_FUS_Controls), colnames(counts2_SC))
counts2_SC_FUS_Controls <- counts2_SC[, common_samples_SC_FUS_Controls]

common_samples_SN_FUS_Controls <- intersect(rownames(metadata_SN_FUS_Controls), colnames(counts2_SN))
counts2_SN_FUS_Controls <- counts2_SN[, common_samples_SN_FUS_Controls]

common_samples_FUS_Controls <- intersect(rownames(metadata_FUS_Controls), colnames(counts2))
counts2_FUS_Controls <- counts2[, common_samples_SC_FUS_Controls]


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


common_samples_SC_ChATneg_FUS_Controls <- intersect(rownames(metadata_SC_ChATneg_FUS_Controls), colnames(counts2_SC_ChATneg))
counts2_SC_ChATneg_FUS_Controls <- counts2_SC_ChATneg[, common_samples_SC_ChATneg_FUS_Controls]

common_samples_SN_ChATneg_FUS_Controls <- intersect(rownames(metadata_SN_ChATneg_FUS_Controls), colnames(counts2_SN_ChATneg))
counts2_SN_ChATneg_FUS_Controls <- counts2_SN_ChATneg[, common_samples_SN_ChATneg_FUS_Controls]








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



setwd("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString_DaCruz_Alejandro/Maria/Analysis/Results/Count_Matrices/Raw/")

write_xlsx(counts_long,"Count_Matrix_Long_Whole_Dataset.xlsx")
write_xlsx(counts_long_SN,"Count_Matrix_Long_Schiatic_Nerve.xlsx")
write_xlsx(counts_long_SC,"Count_Matrix_Long_Spinal_Cord.xlsx")





#Plot the Distributions.

plot <- counts_long %>%
          ggplot(aes(cts, colour = sample)) + 
          geom_freqpoly(binwidth = 1) +
          xlim(0,40)
plot
path <- setwd("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString - DaCruz_Alejandro/Maria/Analysis/Results/Exploratory_Analysis/")
ggsave(filename = "plot.png", plot = plot, device = "png", path = path, width = 12, height = 6)

#Plot the Distributions SN.

plot <- counts_long_SN %>%
  ggplot(aes(cts, colour = sample)) + 
  geom_freqpoly(binwidth = 1) +
  xlim(0,40)
plot
path <- setwd("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString_DaCruz_Alejandro/Maria/Analysis/Results/Exploratory_Analysis/SN")
ggsave(filename = "plot.png", plot = plot, device = "png", path = path, width = 12, height = 6)

#Plot the Distributions SC.

plot <- counts_long_SC %>%
  ggplot(aes(cts, colour = sample)) + 
  geom_freqpoly(binwidth = 1) +
  xlim(0,40)
plot
path <- setwd("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString_DaCruz_Alejandro/Maria/Analysis/Results/Exploratory_Analysis/SC")
ggsave(filename = "plot.png", plot = plot, device = "png", path = path, width = 12, height = 6)


# Plot Distributions by Slide (SN Samples).
plot <- counts_long_SN %>%
  ggplot(aes(cts, colour = `slide name`)) + 
  geom_freqpoly(binwidth = 1) +
  xlim(0,40)
plot
path <- setwd("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString_DaCruz_Alejandro/Maria/Analysis/Results/Exploratory_Analysis/SN")
ggsave(filename = "plot.png", plot = plot, device = "png", path = path, width = 12, height = 6)


# Plot Distributions by Slide (SC Samples).
plot <- counts_long_SC %>%
ggplot(aes(cts, colour = `slide name`)) + 
  geom_freqpoly(binwidth = 1) +
  xlim(0,40)
plot
path <- setwd("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString_DaCruz_Alejandro/Maria/Analysis/Results/Exploratory_Analysis/SN")
ggsave(filename = "plot.png", plot = plot, device = "png", path = path, width = 12, height = 6)

# Plot Distributions by Replicate_ID (SN Samples).
plot <- counts_long_SN %>%
  ggplot(aes(cts, colour = `replicate_id`)) + 
  geom_freqpoly(binwidth = 1) +
  xlim(0,40)
plot
path <- setwd("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString_DaCruz_Alejandro/Maria/Analysis/Results/Exploratory_Analysis/SN")
ggsave(filename = "plot.png", plot = plot, device = "png", path = path, width = 12, height = 6)


# ...............................................................................

# Multiple Linear Regression with Interaction.

# Expression,i = β0,i + β1,i"+",j + β2,i"slide",j + β3,i"+""slide",j
# y = a + b1x1 + b2x2 + b3xix2

# i : gene
# j : sample
# "+" : Chat_pos or Chat_neg
# "slide" : Slide name (ex. CTLR 1,CTRL 2, FUS 3, FUS 4, TDP-43 1, TDP-43 4....)

# Modify the Data.

counts_long_SN$Genes_Sample <- paste(counts_long_SN$Genes, "_", counts_long_SN$sample, sep = "")
counts_long_SN$slide_name <- counts_long_SN$`slide name`
counts_long_SN$slide_name <- str_replace_all(counts_long_SN$slide_name, " ", "_")

# All the Genes Together.

results_simple_lm_sn <- lm(cts ~ slide_name, data = counts_long_SN)
summary(results_simple_lm_sn)

results_simple_lm_rid <- lm(cts ~ replicate_id, data = counts_long_SN)
summary(results_simple_lm_rid)

results_multi_lm <- lm(cts ~ replicate_id + slide_name, data = counts_long_SN)
summary(results_multi_lm)


# Formatting Regression Output.

# 1. With Broom.

library("broom")

tidy_results_simple_lm_sn <-  tidy(results_simple_lm_sn)
tidy_results_simple_lm_rid <-  tidy(results_simple_lm_rid)

# 2. With Stargazer. # Not Working for now.

library("stargazer")

stargazer_results_simple_lm_sn <- stargazer(results_simple_lm_sn, title = "Regression Analysis. Slide Name", out = "tabe1.txt")
stargazer_results_simple_lm_rid <- stargazer(results_simple_lm_rid, header=FALSE)


# 3. Other Ways.

library("jtools")

summ_results_simple_lm_sn <- summ(results_simple_lm_sn)
summ_results_simple_lm_rid <- summ(results_simple_lm_rid)
summ_results_multi_lm <- summ(results_multi_lm)

# 4. Add Coefficient Data to the Initial Dataframe.

# ...... Extract Fitted Values.......
# For the Simple Regression.

fitted_values_rid <- summ_results_simple_lm_rid$model$fitted.values
class(fitted_values_rid)
fv_df <- data.frame(fitted_values_rid)

counts_long_SN_simple_rid_fitted <- cbind(counts_long_SN, fv_df)


fitted_values_sn <- summ_results_simple_lm_sn$model$fitted.values
class(fitted_values_sn)
fv_df <- data.frame(fitted_values_sn)

counts_long_SN_simple_sn_fitted <- cbind(counts_long_SN, fv_df)



# For the Multi Regression.

fitted.values <- summ_results_multi_lm$model$fitted.values
class(fitted.values)
fv_df <- data.frame(fitted.values)

counts_long_SN_fitted <- cbind(counts_long_SN, fv_df)



# Plot Fitted Values.

# For the Single Regression. 

plot <- counts_long_SN_simple_rid_fitted %>%
  ggplot(aes(fitted_values_sn, colour = `replicate_id`)) + 
  geom_freqpoly(binwidth = 1) +
  xlim(0,20)
plot

hist(counts_long_SN_simple_rid_fitted[counts_long_SN_simple_rid_fitted$replicate_id == "FUS_SN_ChaTpos"])

path <- setwd("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString_DaCruz_Alejandro/Maria/Analysis/Results/Exploratory_Analysis/SN")
ggsave(filename = "plot1.png", plot = plot, device = "png", path = path, width = 12, height = 6)


# Plot Distributions by Slide (SN Samples).

plot <- counts_long_SN_fitted %>%
  ggplot(aes(fitted.values, colour = `slide_name`)) + 
  geom_freqpoly(binwidth = 1) +
  xlim(0,20)
plot
path <- setwd("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString_DaCruz_Alejandro/Maria/Analysis/Results/Exploratory_Analysis/SN")
ggsave(filename = "plot2.png", plot = plot, device = "png", path = path, width = 12, height = 6)



# Plot Distributions by Replicate_ID (SN Samples).

plot <- counts_long_SN_fitted %>%
  ggplot(aes(fitted.values, colour = `replicate_id`)) + 
  geom_freqpoly(binwidth = 1) +
  xlim(0,40)
plot
path <- setwd("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString_DaCruz_Alejandro/Maria/Analysis/Results/Exploratory_Analysis/SN")
ggsave(filename = "plot.png", plot = plot, device = "png", path = path, width = 12, height = 6)

# Fitted Values Plot.

ggplot(counts_long_SN, aes(x=predict(results_multi_lm), y = cts)) + 
  geom_point() +
  geom_abline(intercept=0, slope=1) +
  labs(x = ' Predicted Values (Fitted) ', y = ' Actual Values (Initial Counts) ', title='Predicted vs. Actual Values')


# Check for NA in the Data Frame.
sum(is.na(counts_long_SN))
x_nonum <- which(is.na(counts_long_SN$`Genes`))
x_nonum




# ...............................................................................



# Separately per Gene.


Genes <- counts_long_SN$Genes

for(i in 1:Genes){
  Genes <- counts_long_SN$Genes
  results_multi_lm_genes[i] <- lm(cts[i] ~ replicate_id[i], data = counts_long_SN[i])
}
  

# ...............................................................................

# Check the max Counts.

unique(counts_long_SC$cts)
unique(assayData(demoData)[["exprs"]])

Max_cts_SC <- colMax(counts_long_SC)
Max_cts_SC
Max_cts_SN <- colMax(counts_long_SN)
Max_cts_SN

counts_long_SC_over_40_cts <- subset(counts_long_SC, cts > 40)
counts_long_SN_over_40_cts <- subset(counts_long_SN, cts > 40)

counts_long_SC_over_100_cts <- subset(counts_long_SC, cts > 100)
counts_long_SN_over_100_cts <- subset(counts_long_SN, cts > 100)

counts_long_SC_over_500_cts <- subset(counts_long_SC, cts > 500)
counts_long_SN_over_500_cts <- subset(counts_long_SN, cts > 500)

counts_long_SC_over_1000_cts <- subset(counts_long_SC, cts > 1000)
counts_long_SN_over_1000_cts <- subset(counts_long_SN, cts > 1000)




# Save Count Matrices.
setwd("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString_DaCruz_Alejandro/Maria/Analysis/Results/Count_Matrices/Raw/Whole_Dataset")
write_xlsx(counts2, "Count_Matrix.xlsx")


# Save Metadata.
setwd("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString_DaCruz_Alejandro/Maria/Analysis/Results/Metadata")
write_xlsx(metadata, "Metadata.xlsx")







#...................................... Extra.

# Segment QC.


QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (80%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 1,   # Minimum negative control counts (10)
       maxNTCCount = 1000,     # Maximum counts observed in NTC well (1000)
       minNuclei = 20,         # Minimum # of nuclei estimated (100)
       minArea = 1000)         # Minimum segment area (5000)
demoData <-
  setSegmentQCFlags(demoData, 
                    qcCutoffs = QC_params)        

# Collate QC Results
QCResults <- protocolData(demoData)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))





library(ggplot2)

col_by <- "segment"

# Graphical summaries of QC statistics plot function
QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}

QC_histogram(sData(demoData), "Trimmed (%)", col_by, 80)
QC_histogram(sData(demoData), "Stitched (%)", col_by, 80)
QC_histogram(sData(demoData), "Aligned (%)", col_by, 75)
QC_histogram(sData(demoData), "Saturated (%)", col_by, 50) +
  labs(title = "Sequencing Saturation (%)",
       x = "Sequencing Saturation (%)")
QC_histogram(sData(demoData), "area", col_by, 1000, scale_trans = "log10")
QC_histogram(sData(demoData), "nuclei", col_by, 20)
negativeGeoMeans <- 
  esBy(negativeControlSubset(demoData), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 



# calculate the negative geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(demoData), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       })
protocolData(demoData)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData
negCols <- paste0("NegGeoMean_", modules)
pData(demoData)[, negCols] <- sData(demoData)[["NegGeoMean"]]
for(ann in negCols) {
  plt <- QC_histogram(pData(demoData), ann, col_by, 2, scale_trans = "log10")
  print(plt)
}


# detatch neg_geomean columns ahead of aggregateCounts call
pData(demoData) <- pData(demoData)[, !colnames(pData(demoData)) %in% negCols]

# show all NTC values, Freq = # of Segments with a given NTC count:
kable(table(NTC_Count = sData(demoData)$NTC),
      col.names = c("NTC Count", "# of Segments"))


kable(QC_Summary, caption = "QC Summary Table for each Segment")



# Probe QC.


# Generally keep the qcCutoffs parameters unchanged. Set removeLocalOutliers to 
# FALSE if you do not want to remove local outliers
demoData <- setBioProbeQCFlags(demoData, 
                               qcCutoffs = list(minProbeRatio = 0.1,
                                                percentFailGrubbs = 20), 
                               removeLocalOutliers = TRUE)

ProbeQCResults <- fData(demoData)[["QCFlags"]]

# Define QC table for Probe QC
qc_df <- data.frame(Passed = sum(rowSums(ProbeQCResults[, -1]) == 0),
                    Global = sum(ProbeQCResults$GlobalGrubbsOutlier),
                    Local = sum(rowSums(ProbeQCResults[, -2:-1]) > 0
                                & !ProbeQCResults$GlobalGrubbsOutlier))
