library("writexl")
library("DESeq2")
library("apeglm")
library("dplyr")
library("tidyr")
library("reshape2")
library("readxl")


metadata_SN_ChATpos <- metadata_SN_ChATpos[order(match(metadata_SN_ChATpos$sample, colnames(counts2_SN_ChATpos))), ]

dds <- DESeqDataSetFromMatrix(countData = counts2_SN_ChATpos,
                              colData = metadata_SN_ChATpos,
                              design= ~ class)

dds

# Relevel.

dds$class <- relevel(dds$class, "FUS")


# DESeq2 Model.


dds <- DESeq(dds)

assays(dds)[["cooks"]]
par(mar = c(8,5,2,2))
boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2)


# List the Coefficients.

resultsNames(dds)
res <- results(dds, independentFiltering = FALSE)
res <- results(dds)
res
summary(res)


plotDispEsts(dds)
plotMA(res, ylim=c(-2,2))

metadata(res)$alpha
metadata(res)$filterThreshold

plot(metadata(res)$filterNumRej, 
     type = "b", ylab = "number of rejections",
     xlab = "quantiles of filter")
lines(metadata(res)$lo.fit, col = "red")
abline(v = metadata(res)$filterTheta)



plot(res$baseMean+1, -log10(res$pvalue),
     log = "x", xlab = "mean of normalized counts",
     ylab = expression(-log[10](pvalue)),
     ylim = c(0,5),
     #xlim=c(0,20),
     cex = .4, col = rgb(0,0,0,.3))


hist( res$pvalue, breaks=20, col="grey" )



df <- as.data.frame(Results_FUS_Vs_CTRL)
df$Genes <- row.names(df)



setwd("")
write_xlsx(df, "Results_SN_ChATpos_FUS_Vs_CTRL_2.xlsx")


# Results Table and Log fold Change Shrinkage for Visualization and Ranking.



# 0. Region SN Vs SC.

Results_SN_Vs_SC <- results(dds, contrast = c("region", "SN", "SC"))
Results_SN_Vs_SC
Results_SN_Vs_SC_LFC <- lfcShrink(dds, coef = "region_SN_vs_SC", res = Results_SN_Vs_SC , type = "apeglm")
Results_SN_Vs_SC_LFC



# 1. TDP 43 Vs CTRL.

Results_TDP_43_Vs_CTRL <- results(dds, contrast = c("class", "TDP_43", "CTRL"))
Results_TDP_43_Vs_CTRL_LFC <- lfcShrink(dds, coef = "class_TDP_43_vs_CTRL", res = Results_TDP_43_Vs_CTRL , type="apeglm")
Results_TDP_43_Vs_CTRL_LFC



# 2. FUS Vs CTRL.

Results_FUS_Vs_CTRL <- results(dds, contrast = c("class", "FUS", "CTRL"))
Results_FUS_Vs_CTRL_LFC <- lfcShrink(dds, coef = "class_FUS_vs_CTRL", res = Results_FUS_Vs_CTRL , type="apeglm")
Results_FUS_Vs_CTRL_LFC

# 2. TDP 43 Vs FUS.

Results_TDP_43_Vs_FUS <- results(dds, contrast = c("class", "TDP_43", "FUS"))
Results_TDP_43_Vs_FUS_LFC <- lfcShrink(dds, coef = "class_TDP_43_vs_FUS", res = Results_TDP.43_Vs_FUS , type="apeglm")
Results_TDP_43_Vs_FUS_LFC


# 2. ChATpos Vs ChATneg.

Results_ChATpos_Vs_ChATneg <- results(dds, contrast = c("segment", "ChATpos", "ChATneg"))
Results_ChATpos_Vs_ChATneg_LFC <- lfcShrink(dds, coef = "segment_ChATpos_vs_ChATneg", res = Results_ChATpos_Vs_ChATneg , type="apeglm")
Results_ChATpos_Vs_ChATneg_LFC



# We Can Order our Results Table by the Smallest p value.

Results_SN_Vs_SC_LFC_Ordered <- Results_SN_Vs_SC_LFC[order(Results_SN_Vs_SC_LFC$pvalue),]
sum(Results_SN_Vs_SC_LFC$pvalue < 0.05, na.rm=TRUE)

df <- as.data.frame(Results_SN_Vs_SC)
df$Genes <- row.names(df)

df <- as.data.frame(Results_SN_Vs_SC_LFC_Ordered)
df$Genes <- row.names(df)

setwd("")
write_xlsx(df, "Results_SN_Vs_SC_ChATpos_CTRL.xlsx")

#................................................................................................


Results_ChATpos_Vs_ChATneg_Ordered <- Results_ChATpos_Vs_ChATneg[order(Results_ChATpos_Vs_ChATneg$pvalue),]
sum(Results_ChATpos_Vs_ChATneg_Ordered$pvalue < 0.05, na.rm=TRUE)

df <- as.data.frame(Results_ChATpos_Vs_ChATneg_Ordered)
df$Genes <- row.names(df)

setwd("")
write_xlsx(df, "Results_ChATpos_Vs_ChATneg_Ordered.xlsx")


#................................................................................................


Results_TDP_43_Vs_CTRL_Ordered <- Results_TDP_43_Vs_CTRL[order(Results_TDP_43_Vs_CTRL$pvalue),]
sum(Results_TDP_43_Vs_CTRL_Ordered$pvalue < 0.05, na.rm=TRUE)

df <- as.data.frame(Results_TDP_43_Vs_CTRL_Ordered)
df$Genes <- row.names(df)

setwd("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString_DaCruz_Alejandro/Maria/Analysis/Results/DESeq2_Analysis/SC/Excel_Files")
write_xlsx(df, "Results_TDP_43_Vs_CTRL_Ordered.xlsx")


#................................................................................................

Results_FUS_Vs_CTRL_Ordered <- Results_FUS_Vs_CTRL[order(Results_FUS_Vs_CTRL$pvalue),]
sum(Results_FUS_Vs_CTRL_Ordered$pvalue < 0.05, na.rm=TRUE)

df <- as.data.frame(Results_FUS_Vs_CTRL_Ordered)
df$Genes <- row.names(df)

setwd("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString_DaCruz_Alejandro/Maria/Analysis/Results/DESeq2_Analysis/SC/Excel_Files")
write_xlsx(df, "Results_FUS_Vs_CTRL_Ordered.xlsx")

# ...........................................................................................

Results_TDP_43_Vs_FUS_Ordered <- Results_TDP_43_Vs_FUS[order(Results_TDP_43_Vs_FUS$pvalue),]
sum(Results_TDP_43_Vs_FUS_Ordered$pvalue < 0.05, na.rm=TRUE)

df <- as.data.frame(Results_TDP_43_Vs_FUS_Ordered)
df$Genes <- row.names(df)

setwd("C:/Users/u0140675/Desktop/Μαρία/PhD/Projects/NanoString_DaCruz_Alejandro/Maria/Analysis/Results/DESeq2_Analysis/SC/Excel_Files")
write_xlsx(df, "Results_TDP_43_Vs_FUS_Ordered.xlsx")


