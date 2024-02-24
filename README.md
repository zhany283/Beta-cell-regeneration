# Beta-cell-regeneration
library(Seurat)
library(decontX)
library(ggplot2)
library(dplyr)
library(Matrix)
library(SeuratObject)

setwd("C:\\Users\\ethan\\Desktop\\solomon_scRNA-Seq\\")


FoxO1_KO_1 <- Read10X(data.dir = "FoxO1_KO_1/filtered_feature_bc_matrix")
FoxO1_KO_2 <- Read10X(data.dir = "FoxO1_KO_2/filtered_feature_bc_matrix")
HCD_1 <- Read10X(data.dir = "HCD_1/filtered_feature_bc_matrix")
HCD_2 <- Read10X(data.dir = "HCD_2/filtered_feature_bc_matrix")
WT_1 <- Read10X(data.dir = "WT_1/filtered_feature_bc_matrix")
WT_2 <- Read10X(data.dir = "WT_2/filtered_feature_bc_matrix")



FoxO1_KO_1_seurat <- CreateSeuratObject(counts = FoxO1_KO_1, project = "FoxO1_KO_1")
FoxO1_KO_2_seurat <- CreateSeuratObject(counts = FoxO1_KO_2, project = "FoxO1_KO_2")
HCD_1_seurat <- CreateSeuratObject(counts = HCD_1, project = "HCD_1")
HCD_2_seurat <- CreateSeuratObject(counts = HCD_2, project = "HCD_2")
WT_1_seurat <- CreateSeuratObject(counts = WT_1, project = "WT_1")
WT_2_seurat <- CreateSeuratObject(counts = WT_2, project = "WT_2")


counts_FoxO1_KO_1 <- FoxO1_KO_1_seurat@assays$RNA@layers$counts
counts_FoxO1_KO_2 <- FoxO1_KO_2_seurat@assays$RNA@layers$counts
counts_HCD_1 <- HCD_1_seurat@assays$RNA@layers$counts
counts_HCD_2 <- HCD_2_seurat@assays$RNA@layers$counts
counts_WT_1 <- WT_1_seurat@assays$RNA@layers$counts
counts_WT_2 <- WT_2_seurat@assays$RNA@layers$counts


decontX_results_FoxO1_KO_1 <- decontX(counts_FoxO1_KO_1)
decontX_results_FoxO1_KO_2 <- decontX(counts_FoxO1_KO_2)
decontX_results_HCD_1 <- decontX(counts_HCD_1)
decontX_results_HCD_2 <- decontX(counts_HCD_2)
decontX_results_WT_1 <- decontX(counts_WT_1)
decontX_results_WT_2 <- decontX(counts_WT_2)

FoxO1_KO_1_seurat$Contamination <- decontX_results_FoxO1_KO_1$contamination
FoxO1_KO_2_seurat$Contamination <- decontX_results_FoxO1_KO_2$contamination
HCD_1_seurat$Contamination <- decontX_results_HCD_1$contamination
HCD_2_seurat$Contamination <- decontX_results_HCD_2$contamination
WT_1_seurat$Contamination <- decontX_results_WT_1$contamination
WT_2_seurat$Contamination <- decontX_results_WT_2$contamination


low_con_FoxO1_KO_1_seurat <- FoxO1_KO_1_seurat[, FoxO1_KO_1_seurat$Contamination < 0.2]
low_con_FoxO1_KO_2_seurat <- FoxO1_KO_2_seurat[, FoxO1_KO_2_seurat$Contamination < 0.2]
low_con_HCD_1_seurat <- HCD_1_seurat[, HCD_1_seurat$Contamination < 0.2]
low_con_HCD_2_seurat <- HCD_2_seurat[, HCD_2_seurat$Contamination < 0.2]
low_con_WT_1_seurat <- WT_1_seurat[, WT_1_seurat$Contamination < 0.2]
low_con_WT_2_seurat <- WT_2_seurat[, WT_2_seurat$Contamination < 0.2]

seuratObject <- merge(low_con_FoxO1_KO_1_seurat, y = list(low_con_FoxO1_KO_2_seurat,
                                                          low_con_HCD_1_seurat,
                                                          low_con_HCD_2_seurat,
                                                          low_con_WT_1_seurat,                                                             
                                                          low_con_WT_2_seurat),
                      add.cell.ids = c("FoxO1_KO_1", "FoxO1_KO_2", "HCD_1", "HCD_2", "WT_1", "WT_2"))
#save RDS
saveRDS(seuratObject, file = "C:\\Users\\ethan\\Desktop\\solomon_scRNA-Seq\\WTvsHCDvsFoxo1_02222024.rds")
#load back
seuratObject <- readRDS(file = "C:\\Users\\ethan\\Desktop\\solomon_scRNA-Seq\\WTvsHCDvsFoxo1_02222024.rds")
# The pattern might need to be adjusted depending on your organism and gene naming conventions
mitochondrial_genes <- grep("^MT-", rownames(seuratObject), value = TRUE)

# Calculate the percentage of mitochondrial gene expression
seuratObject[["percent.mt"]] <- PercentageFeatureSet(seuratObject, features = mitochondrial_genes)

#QC
seuratObject <- subset(seuratObject, subset = nFeature_RNA > 500 & nFeature_RNA < 2500 & percent.mt < 10)

#Normalize the data
seuratObject <- NormalizeData(seuratObject, normalization.method = "LogNormalize", scale.factor = 10000)
#Find Variable Features
seuratObject <- FindVariableFeatures(seuratObject, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seuratObject), 10)
plot1 <- VariableFeaturePlot(seuratObject)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
![VariableFeaturePlot](https://github.com/zhany283/Beta-cell-regeneration/assets/130387837/bb4c88a5-c6da-461a-9360-ec0d3d82a91c)

