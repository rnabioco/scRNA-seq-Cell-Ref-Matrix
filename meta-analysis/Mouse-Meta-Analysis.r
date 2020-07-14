library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(here)

# figure out project root
proj_dir <- here()
mouseAtlas <- readRDS(file.path(proj_dir, "atlas", "musMusculus", "mouseAtlas.rds"))

mouseMetaAnalysis <- CreateSeuratObject(counts = mouseAtlas, project = "Mouse-Meta-Analysis", min.cells = 0, min.features = 0)
gc()

#Normalize Data
mouseMetaAnalysis <- NormalizeData(mouseMetaAnalysis, normalization.method = "LogNormalize", scale.factor = 10000)

#Preprocessing workflow
#mouseMetaAnalysis@assays$RNA@data <- mouseMetaAnalysis@assays$RNA@counts
mouseMetaAnalysis[["percent.mt"]] <- PercentageFeatureSet(mouseMetaAnalysis, pattern = "^mt-")
head(mouseMetaAnalysis@meta.data, 5)
VlnPlot(mouseMetaAnalysis, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(mouseMetaAnalysis, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(mouseMetaAnalysis, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Find Variable Features for PCA
mouseMetaAnalysis <- FindVariableFeatures(mouseMetaAnalysis, selection.method = "mean.var.plot", nfeatures = 2000)
top10 <- head(VariableFeatures(mouseMetaAnalysis), 10)
plot1 <- VariableFeaturePlot(mouseMetaAnalysis)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Linear dimension reduction/Run PCA
all.genes <- rownames(mouseMetaAnalysis)
mouseMetaAnalysis <- ScaleData(mouseMetaAnalysis, features = all.genes)
#mouseMetaAnalysis <- RunPCA(mouseMetaAnalysis, features = VariableFeatures(object = mouseMetaAnalysis), npcs = 49)
mouseMetaAnalysis <- RunPCA(mouseMetaAnalysis,
                            features = VariableFeatures(object = mouseMetaAnalysis),
                            npcs = ncol(mouseMetaAnalysis) - 1)
print(mouseMetaAnalysis[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(mouseMetaAnalysis, dims = 1:2, reduction = "pca")
DimPlot(mouseMetaAnalysis, reduction = "pca")
DimHeatmap(mouseMetaAnalysis, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(mouseMetaAnalysis, dims = 1:15, cells = 500, balanced = TRUE)

#Determine dimensionality
mouseMetaAnalysis <- JackStraw(mouseMetaAnalysis, num.replicate = 100)
mouseMetaAnalysis <- ScoreJackStraw(mouseMetaAnalysis, dims = 1:20)
JackStrawPlot(mouseMetaAnalysis, dims = 1:15)
ElbowPlot(mouseMetaAnalysis)

#Clustering
mouseMetaAnalysis <- FindNeighbors(mouseMetaAnalysis, dims = 1:10)
mouseMetaAnalysis <- FindClusters(mouseMetaAnalysis, resolution = 10.0)
head(Idents(mouseMetaAnalysis), 5)

#Create unannotated UMAP
mouseMetaAnalysis <- RunUMAP(mouseMetaAnalysis, dims = 1:10)
DimPlot(mouseMetaAnalysis, reduction = "umap")

#Differentially expressed features
cluster1.markers <- FindMarkers(mouseMetaAnalysis, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
cluster5.markers <- FindMarkers(mouseMetaAnalysis, ident.1 = 5, ident.2 = c(0,3), min.pct = 0.25)
head(cluster5.markers, n = 5)
mouseMetaAnalysis.markers <- FindAllMarkers(mouseMetaAnalysis, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
mouseMetaAnalysis.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
cluster1.markers <- FindMarkers(mouseMetaAnalysis, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
top10 <- mouseMetaAnalysis.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC) #Create top 10 markers for each cluster
DoHeatmap(pbmc, features = top10$gene) + NoLegend() #Create heat map of top 10 markers

#Assign cell types/Annotate UMAP
new.cluster.ids <- colnames(mouseAtlas)
names(new.cluster.ids) <- levels(mouseMetaAnalysis)
mouseMetaAnalysis <- RenameIdents(mouseMetaAnalysis, new.cluster.ids)
AnnotatedUMAP <- DimPlot(mouseMetaAnalysis, reduction = "umap", label = TRUE, pt.size = 0.5) 
HoverLocator(plot = AnnotatedUMAP, information = FetchData(mouseMetaAnalysis, vars = c("seurat_clusters")))
mouseMetaAnalysis@meta.data$study <- str_remove(rownames(mouseMetaAnalysis@meta.data), ".+\\(") %>% str_remove("\\)")
DimPlot(mouseMetaAnalysis, reduction = "umap", group.by = "study")
