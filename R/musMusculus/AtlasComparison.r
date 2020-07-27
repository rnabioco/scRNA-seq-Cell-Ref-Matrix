library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(digest)
library(here)
library(Matrix)

proj_dir <- here()
proj_dir <- file.path("Reference-Matrix-Generation", "data", "atlasComparison", "GSE141259")

mat_Lung <- Read10X(data.dir = proj_dir, gene.column = 1)
SeuratLung <- CreateSeuratObject(counts = mat_Lung, project = "mat_Lung", min.cells = 3, min.features = 200)
SeuratLung

SeuratLung[["percent.mt"]] <- PercentageFeatureSet(SeuratLung, pattern = "^mt-")
# Visualize QC metrics as a violin plot
VlnPlot(SeuratLung, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(SeuratLung, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SeuratLung, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
