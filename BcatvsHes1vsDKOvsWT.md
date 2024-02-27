library(Seurat)
library(decontX)
library(ggplot2)
library(dplyr)
library(Matrix)
library(SeuratObject)

setwd("C:\\Users\\ethan\\Desktop\\solomon_scRNA-Seq\\")

bCatenin_KO_1 <- Read10X(data.dir = "bCatenin_KO_1/filtered_feature_bc_matrix")
bCatenin_KO_2 <- Read10X(data.dir = "bCatenin_KO_2/filtered_feature_bc_matrix")
DKO_1 <- Read10X(data.dir = "DKO_1/filtered_feature_bc_matrix")
DKO_2 <- Read10X(data.dir = "DKO_2/filtered_feature_bc_matrix")
Hes1_KO_1 <- Read10X(data.dir = "Hes1_KO1/filtered_feature_bc_matrix")
Hes1_KO_2 <- Read10X(data.dir = "Hes1_KO2/filtered_feature_bc_matrix")
WT_1 <- Read10X(data.dir = "WT_1/filtered_feature_bc_matrix")
WT_2 <- Read10X(data.dir = "WT_2/filtered_feature_bc_matrix")


bCatenin_KO_1_seurat <- CreateSeuratObject(counts = bCatenin_KO_1, project = "bCatenin_KO_1")
bCatenin_KO_2_seurat <- CreateSeuratObject(counts = bCatenin_KO_2, project = "bCatenin_KO_2")
DKO_1_seurat <- CreateSeuratObject(counts = DKO_1, project = "DKO_1")
DKO_2_seurat <- CreateSeuratObject(counts = DKO_2, project = "DKO_2")
Hes1_KO_1_seurat <- CreateSeuratObject(counts = Hes1_KO_1, project = "Hes1_KO1")
Hes1_KO_2_seurat <- CreateSeuratObject(counts = Hes1_KO_2, project = "Hes1_KO2")
WT_1_seurat <- CreateSeuratObject(counts = WT_1, project = "WT_1")
WT_2_seurat <- CreateSeuratObject(counts = WT_2, project = "WT_2")

counts_bCatenin_KO_1 <- bCatenin_KO_1_seurat@assays$RNA@layers$counts
counts_bCatenin_KO_2 <- bCatenin_KO_2_seurat@assays$RNA@layers$counts
counts_DKO_1 <- DKO_1_seurat@assays$RNA@layers$counts
counts_DKO_2 <- DKO_2_seurat@assays$RNA@layers$counts
counts_Hes1_KO1 <- Hes1_KO_1_seurat@assays$RNA@layers$counts
counts_Hes1_KO2 <- Hes1_KO_2_seurat@assays$RNA@layers$counts
counts_WT_1 <- WT_1_seurat@assays$RNA@layers$counts
counts_WT_2 <- WT_2_seurat@assays$RNA@layers$counts

decontX_results_bCatenin_KO_1 <- decontX(counts_bCatenin_KO_1) 
decontX_results_bCatenin_KO_2 <- decontX(counts_bCatenin_KO_2)
decontX_results_DKO_1 <- decontX(counts_DKO_1)
decontX_results_DKO_2 <- decontX(counts_DKO_2)
decontX_results_Hes1_KO1 <- decontX(counts_Hes1_KO1)
decontX_results_Hes1_KO2 <- decontX(counts_Hes1_KO2)
decontX_results_WT_1 <- decontX(counts_WT_1)
decontX_results_WT_2 <- decontX(counts_WT_2)

bCatenin_KO_1_seurat$Contamination <- decontX_results_bCatenin_KO_1$contamination
bCatenin_KO_2_seurat$Contamination <- decontX_results_bCatenin_KO_2$contamination
DKO_1_seurat$Contamination <- decontX_results_DKO_1$contamination
DKO_2_seurat$Contamination <- decontX_results_DKO_2$contamination
Hes1_KO_1_seurat$Contamination <- decontX_results_Hes1_KO1$contamination
Hes1_KO_2_seurat$Contamination <- decontX_results_Hes1_KO2$contamination
WT_1_seurat$Contamination <- decontX_results_WT_1$contamination
WT_2_seurat$Contamination <- decontX_results_WT_2$contamination


low_con_bCatenin_KO_1_seurat= bCatenin_KO_1_seurat[,bCatenin_KO_1_seurat$Contamination < 0.2]
low_con_bCatenin_KO_2_seurat <- bCatenin_KO_2_seurat[, bCatenin_KO_2_seurat$Contamination < 0.2]
low_con_DKO_1_seurat <- DKO_1_seurat[, DKO_1_seurat$Contamination < 0.2]
low_con_DKO_2_seurat <- DKO_2_seurat[, DKO_2_seurat$Contamination < 0.2]
low_con_Hes1_KO_1_seurat <- Hes1_KO_1_seurat[, Hes1_KO_1_seurat$Contamination < 0.2]
low_con_Hes1_KO_2_seurat <- Hes1_KO_2_seurat[, Hes1_KO_2_seurat$Contamination < 0.2]
low_con_WT_1_seurat <- WT_1_seurat[, WT_1_seurat$Contamination < 0.2]
low_con_WT_2_seurat <- WT_2_seurat[, WT_2_seurat$Contamination < 0.2]

seuratObject <- merge(low_con_bCatenin_KO_1_seurat, y = list(low_con_bCatenin_KO_2_seurat,
                                                             low_con_DKO_1_seurat,
                                                             low_con_DKO_2_seurat,
                                                             low_con_Hes1_KO_1_seurat,
                                                             low_con_Hes1_KO_2_seurat,
                                                             low_con_WT_1_seurat,                                                             
                                                             low_con_WT_2_seurat),
                      add.cell.ids = c("bCatenin_KO_1", "bCatenin_KO_2", "DKO_1", "DKO_2", "Hes1_KO1", "Hes1_KO2", "WT_1", "WT_2"))
#save RDS
saveRDS(seuratObject, file = "C:\\Users\\ethan\\Desktop\\solomon_scRNA-Seq\\seurat_object_bcat.rds")

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
![image](https://github.com/zhany283/Beta-cell-regeneration/assets/130387837/38ae7ded-c968-4319-ac5a-0854b8a78b16)

#scale the data
seuratObject <- ScaleData(seuratObject, features = VariableFeatures(object = seuratObject))

#Linear dimensional reduction
seuratObject <- RunPCA(seuratObject, features = VariableFeatures(object = seuratObject))

DimHeatmap(seuratObject, dims = 1:10, cells = 500, balanced = TRUE)
![image](https://github.com/zhany283/Beta-cell-regeneration/assets/130387837/969fc1ed-26b3-429b-818e-bc1dca78e802)


#cluster
seuratObject <- FindNeighbors(seuratObject, dims = 1:10)
seuratObject <- FindClusters(seuratObject, resolution = 0.5)

#saveRDS(seuratObject, file = "C:\\Users\\ethan\\Desktop\\solomon_scRNA-Seq\\seuratObject_02222024.rds")

seuratObject <- RunTSNE(seuratObject, dims = 1:10)
DimPlot(seuratObject, reduction = "tsne")

seuratObject <- RunUMAP(seuratObject, dims = 1:10)
DimPlot(seuratObject, reduction = "umap")
![image](https://github.com/zhany283/Beta-cell-regeneration/assets/130387837/abba4e72-b3ad-44a8-b702-ecfc754d327c)


VlnPlot(seuratObject, features = c("tdTomato-all","Krt19","Sox9","Spp1","Muc1","Ins1","Ucn3","Pdx1","Mnx1","Mafa","Nkx6-1","Sst","Hhex","Gcg","Ppy","Hes1","Ins2","Nkx2-2","Pax6"), pt.size = 0.1)
![image](https://github.com/zhany283/Beta-cell-regeneration/assets/130387837/59427775-f199-4bb3-8f3e-31801a5b6807)
VlnPlot(seuratObject, features = c("Cd14","Cd3e","Cd19","Ghrl","Cd14","Mrc1","Vim","Acta2","Fap","Cldn4","Prrx1","Ctrb1","Pecam1","Epcam","Cpa1"), pt.size = 0.1)
![image](https://github.com/zhany283/Beta-cell-regeneration/assets/130387837/358821b0-fb8d-481e-9142-5ba5402034e6)

