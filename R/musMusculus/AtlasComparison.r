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

SeuratLung <- NormalizeData(SeuratLung)
SeuratLung <- FindVariableFeatures(SeuratLung, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SeuratLung), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(SeuratLung)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(SeuratLung)
SeuratLung <- ScaleData(SeuratLung, features = all.genes)
SeuratLung <- RunPCA(SeuratLung, features = VariableFeatures(object = SeuratLung))
print(SeuratLung[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(SeuratLung, dims = 1:2, reduction = "pca")
DimPlot(SeuratLung, reduction = "pca")
DimHeatmap(SeuratLung, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(SeuratLung, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(SeuratLung)
SeuratLung <- FindNeighbors(SeuratLung, dims = 1:10)
SeuratLung <- FindClusters(SeuratLung, resolution = 0.5)
head(Idents(SeuratLung), 5)

SeuratLung <- RunUMAP(SeuratLung, dims = 1:10)
DimPlot(SeuratLung, reduction = "umap")

cluster1.markers <- FindMarkers(SeuratLung, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
cluster5.markers <- FindMarkers(SeuratLung, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
SeuratLung.markers <- FindAllMarkers(SeuratLung, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SeuratLung.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
cluster1.markers <- FindMarkers(SeuratLung, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
top10 <- SeuratLung.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(SeuratLung, features = top10$gene) + NoLegend()
saveRDS(SeuratLung, file = "Lung.rds")
