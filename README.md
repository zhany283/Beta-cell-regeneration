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
 #scale the data
seuratObject <- ScaleData(seuratObject, features = VariableFeatures(object = seuratObject))
Centering and scaling data matrix
  |=============================================================================================================| 100%
#Linear dimensional reduction
seuratObject <- RunPCA(seuratObject, features = VariableFeatures(object = seuratObject))
PC_ 1 
Positive:  Tmsb4x, Tpt1, Rpl18a, Rps11, Rps16, Rpl23, Fau, Ppia, Rps29, Rps4x 
	   Rps12, Rps9, Rps14, Rplp1, Rps8, Rps23, Eef1a1, Rps5, Actb, Rps24 
	   Fth1, Rps20, Rps19, Rps7, Rpl9, Rpl7, Rpl28, Rpl21, Rpl35a, Rps3 
Negative:  Phactr1, Hs6st3, Sgcz, Car10, Cntnap2, Dlgap1, Pde10a, Rgs7, Fhit, Ins1 
	   Kcnb2, Ins2, Esrrg, Fgf12, Frmd5, Trpm3, Iapp, Gm26917, Ctnna2, Pex5l 
	   Kcnip4, 5330434G04Rik, Galnt14, Npas3, Ptprt, Pcdh15, G6pc2, Gm34759, Sdk1, Tnr 
PC_ 2 
Positive:  Laptm5, Tyrobp, Ctss, H2-Aa, Rgs1, H2-Eb1, H2-Ab1, Cd74, Ptprc, Cybb 
	   Fyb, Fcer1g, Cd83, Lyz2, Mir142hg, Il1b, Cd52, Arhgap15, H2-DMa, C1qc 
	   Ly86, Coro1a, Aif1, Cd86, C1qb, Mpeg1, Cd300c2, Ms4a6c, Ms4a7, Bcl2a1b 
Negative:  Flt1, Plvap, Emcn, Plpp3, Ptprb, Egfl7, Kdr, Ly6c1, Eng, Adgrl4 
	   Plpp1, Esam, Sparc, Adgrf5, Col4a1, Cyyr1, Ramp2, Igfbp7, Clec14a, Slc9a3r2 
	   Pecam1, Podxl, Cdh5, Col4a2, Mgll, Mmrn2, Timp3, Cavin3, Epas1, Cd200 
PC_ 3 
Positive:  Srgn, Mef2c, Fxyd5, Msn, AW112010, Flt1, Plvap, Emcn, Kdr, Ptprb 
	   Ehd4, Plekho1, Gngt2, Adgrl4, Eng, Pecam1, Ly6c1, Egfl7, Arhgap31, Nrp1 
	   Cyyr1, Fli1, Entpd1, Gimap6, Kctd12, Tnfaip2, Inpp5d, Plpp1, Adgrf5, Ctla2a 
Negative:  Spp1, Krt18, Cldn3, Fxyd3, Clu, Krt8, Epcam, Atp1b1, Krt19, Tstd1 
	   tdTomato-all, Wfdc2, Mgst1, Cldn7, Tspan8, Bicc1, Dbi, Gsta3, Naaladl2, Pdzk1ip1 
	   Ambp, Cp, Prr15l, Pkhd1, Sdc4, Tm4sf4, Ehf, Wwc1, Smim6, Sox9 
PC_ 4 
Positive:  Clu, Spp1, Cldn3, Cd24a, Atp1b1, Krt8, Krt18, Fxyd3, Epcam, Krt19 
	   Anxa3, Tstd1, Wfdc2, Pdgfd, Hhex, Tspan8, Cldn7, Rgcc, Ly6e, Kdr 
	   tdTomato-all, Flt1, Hbegf, Emcn, Srgn, Plscr2, Rab27b, Tm4sf4, Smim6, Prr15l 
Negative:  Col3a1, Col1a2, Dcn, Cygb, Rarres2, Col6a1, Serping1, Fbln1, C1s1, Col1a1 
	   Col6a2, Cfh, Crispld2, Serpinf1, Igfbp6, Lum, Pcolce, Olfml3, Gdf10, Itgbl1 
	   Pdgfra, Prrx1, Abca8a, Mfap4, Mmp2, Lama2, Gsn, Hsd11b1, Nbl1, Htra3 
PC_ 5 
Positive:  C1qc, C1qb, Ms4a7, C1qa, Csf1r, Aif1, Cybb, Ly86, Cd300c2, Lyz2 
	   Adgre1, Ctss, Cx3cr1, Mpeg1, Gatm, Apoe, Trem2, H2-Aa, Fcgr3, H2-Eb1 
	   C3ar1, Ccl3, H2-Ab1, Lgmn, Lair1, Cd74, Slamf9, Lst1, Ccl4, Cd68 
Negative:  Itk, Cd3g, Il7r, Cxcr6, Trbc2, Ptprcap, Ikzf3, Ltb, Icos, Vps37b 
	   S100a4, Cd3d, Gimap3, Cd3e, Hcst, Ptpn22, Itgb7, Bcl11b, Trbc1, Sept1 
	   Camk4, Rac2, Stat4, Lat, Il18r1, Cd226, P2ry10, Samsn1, D16Ertd472e, Podnl1 

DimHeatmap(seuratObject, dims = 1:10, cells = 500, balanced = TRUE)
![image](https://github.com/zhany283/Beta-cell-regeneration/assets/130387837/2fdbc7f8-b61a-423a-8b6c-65c5791339a3)
#cluster
seuratObject <- FindNeighbors(seuratObject, dims = 1:10)
seuratObject <- FindClusters(seuratObject, resolution = 0.5)
#saveRDS(seuratObject, file = "C:\\Users\\ethan\\Desktop\\solomon_scRNA-Seq\\seuratObject_02222024.rds")
seuratObject <- RunUMAP(seuratObject, dims = 1:10)
DimPlot(seuratObject, reduction = "umap", label = TRUE)
![image](https://github.com/zhany283/Beta-cell-regeneration/assets/130387837/fe312c3c-ad1a-4851-9a3e-ae18e5a842af)
![image](https://github.com/zhany283/Beta-cell-regeneration/assets/130387837/febfff2e-6cd6-4f3e-86b5-ec5e5aa9bbb2)

VlnPlot(seuratObject, features = c("tdTomato-all","Krt19","Sox9","Spp1","Muc1","Ins1","Ucn3","Pdx1","Mnx1","Mafa","Nkx6-1","Sst","Hhex","Gcg","Ppy","Hes1","Ins2","Nkx2-2","Pax6"), pt.size = 0.1)
![image](https://github.com/zhany283/Beta-cell-regeneration/assets/130387837/2df7ad27-e069-4682-bd9c-bee7742920f2)
#the plot below need to be revised -- 
VlnPlot(seuratObject, features = c("tdTomato-all","Krt19","Sox9","Spp1","Muc1","Ins1","Ucn3","Pdx1","Mnx1","Mafa","Nkx6-1","Sst","Hhex","Gcg","Ppy","Hes1","Ins2","Nkx2-2","Pax6"), group.by = "orig.ident", pt.size = 0.1)


