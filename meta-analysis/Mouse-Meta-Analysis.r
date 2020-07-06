library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)

mouseAtlas <- readRDS(file.path("~/Reference-Matrix-Generation/atlas/musMusculus/MouseAtlas.rds"))

mouseMetaAnalysis <- CreateSeuratObject(counts = mouseAtlas, project = "Mouse-Meta-Analysis", min.cells = 3, min.features = 200)
mouseMetaAnalysis
gc()


#Preprocessing workflow
mouseMetaAnalysis@assays$RNA@data <- mouseMetaAnalysis@assays$RNA@counts
mouseMetaAnalysis[["percent.mt"]] <- PercentageFeatureSet(mouseMetaAnalysis, pattern = "^MT")
head(mouseMetaAnalysis@meta.data, 5)
VlnPlot(mouseMetaAnalysis, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(mouseMetaAnalysis, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mouseMetaAnalysis, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2