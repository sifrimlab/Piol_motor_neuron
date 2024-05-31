library("readxl")
library("ggplot2")
library("ggplotify")
library("ggrepel")
library("ggforce")
library("ggalt")
library("dplyr")

setwd("C:/Users/u0140675/OneDrive - KU Leuven/DESeq2_Analysis_Nanostring/No_Shrunken/SN/Excel_Files/With_Interaction_Terms/ChATpos/")

SC_ChAT <- read_excel("Results_SN_FUS_Vs_CTRL_ChATpos_Ordered.xlsx")
SC_ChAT <- SC_ChAT[,c("Genes", "log2FoldChange", "pvalue")]
sum(is.na(SC_ChAT$pvalue))
SC_ChAT <- SC_ChAT[complete.cases(SC_ChAT),]
sum(is.na(SC_ChAT$pvalue))


setwd("C:/Users/u0140675/OneDrive - KU Leuven/DESeq2_Analysis_Nanostring/Plots/Volcano_Plots/Without_Shrinkage/Color_Version_2/")

vlc_plot <- ggplot(SC_ChAT, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = ifelse(pvalue > 0.05 | abs(log2FoldChange) < 1, "grey",
                                ifelse(log2FoldChange > 0, "plum", "purple")))) +
  scale_color_manual(values = c("grey", "plum", "purple"), 
                     name = "Significance",
                     labels = c("Non-significant", "FUS", "CTRL")) +
  theme_bw() +
  labs(x = "Log2 Fold Change", 
       y = "-Log10(p-value)", 
       title = "Volcano plot of Differentially Expressed Genes",
       subtitle = "Sciatic Nerve - ChATpos - FUS Vs CTRL",
       color = "Significant") +
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "blue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")

vlc_plot

ggsave(vlc_plot, filename = "Volcano_Plot_SN_ChATpos_FUS_Vs_CTRL.png", width = 8, height = 6)
