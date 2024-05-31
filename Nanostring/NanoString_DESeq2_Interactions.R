library("writexl")
library("DESeq2")
library("apeglm")


dds <- DESeqDataSetFromMatrix(countData = counts2_SN,
                              colData = metadata_SN,
                              design= ~ segment + class + segment:class)

dds


# DESeq2 Model.

dds <- DESeq(dds)

# List the Coefficients.

resultsNames(dds)
res <- results(dds)
res

kable( colData(dds))


# Get the Model Matrix

mod_mat <- model.matrix(design(dds), colData(dds))

# Define Coefficient Vectors for Each Condition.

CTRL_ChATpos <- colMeans(mod_mat[dds$class == "CTRL" & dds$segment == "ChATpos", ])
CTRL_ChATneg <- colMeans(mod_mat[dds$class == "CTRL" & dds$segment == "ChATneg", ])

FUS_ChATpos <- colMeans(mod_mat[dds$class == "FUS" & dds$segment == "ChATpos", ])
FUS_ChATneg <- colMeans(mod_mat[dds$class == "FUS" & dds$segment == "ChATneg", ])



TDP.43_ChATpos <- colMeans(mod_mat[dds$class == "TDP.43" & dds$segment == "ChATpos", ])
TDP.43_ChATneg <- colMeans(mod_mat[dds$class == "TDP.43" & dds$segment == "ChATneg", ])





