library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)

humanAtlas <- readRDS(file.path("~/Reference-Matrix-Generation/atlas/homoSapiens/HumanAtlas.rds"))

humanMetaAnalysis <- CreateSeuratObject(counts = humanAtlas, project = "Human-Meta-Analysis", min.cells = 3, min.features = 200)
humanMetaAnalysis
#rm(humanAtlas)
gc()
humanMetaAnalysis@assays$RNA@data <- humanMetaAnalysis@assays$RNA@counts
humanMetaAnalysis[["percent.mt"]] <- PercentageFeatureSet(humanMetaAnalysis, pattern = "^MT")
head(humanMetaAnalysis@meta.data, 5)
VlnPlot(humanMetaAnalysis, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(humanMetaAnalysis, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(humanMetaAnalysis, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
humanMetaAnalysis <- subset(humanMetaAnalysis, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

humanMetaAnalysis <- NormalizeData(humanMetaAnalysis, normalization.method = "LogNormalize", scale.factor = 10000)

humanMetaAnalysis <- FindVariableFeatures(humanMetaAnalysis, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(humanMetaAnalysis), 10)
plot1 <- VariableFeaturePlot(humanMetaAnalysis)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Linear dimension reduction
all.genes <- rownames(humanMetaAnalysis)
humanMetaAnalysis <- ScaleData(humanMetaAnalysis, features = all.genes)
humanMetaAnalysis <- RunPCA(humanMetaAnalysis, features = VariableFeatures(object = humanMetaAnalysis))
print(humanMetaAnalysis[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(humanMetaAnalysis, dims = 1:2, reduction = "pca")
DimPlot(humanMetaAnalysis, reduction = "pca")
DimHeatMap(humanMetaAnalysis, dims = 1, cells = 500, balanced = TRUE)
DimHeatMap(humanMetaAnalysis, dims = 1:15, cells = 500, balanced = TRUE)

#Dimensionality
ElbowPlot(humanMetaAnalysis)

#Cluster cells
humanMetaAnalysis <- FindNeighbors(humanMetaAnalysis, dims = 1:10)
humanMetaAnalysis <- FindClusters(humanMetaAnalysis, resolution = 0.5)
head(Idents(humanMetaAnalysis), 5)

#Non-linear dimensional reduction (UMAP/tSNE)
humanMetaAnalysis <- RunUMAP(humanMetaAnalysis, dims = 1:10)
DimPlot(humanMetaAnalysis, reduction = "umap")
