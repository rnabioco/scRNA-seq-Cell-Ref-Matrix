library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)

D0_FACSatlas <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE143nnn/GSE143435/suppl/GSE143435_DeMicheli_D0_FACSatlas_normalizeddata.txt.gz")
D0_FACSatlas <- D0_FACSatlas %>%
  as.data.frame() %>%
  column_to_rownames('X1') %>%
  as.matrix() %>%
  t()
D0_FACSatlas[1:5, 1:5]

D0_FACSatlasMetadata <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE143nnn/GSE143435/suppl/GSE143435_DeMicheli_D0_FACSatlas_metadata.txt.gz")
sum(colnames(D0_FACSatlas) %in% D0_FACSatlasMetadata$X1)
ncol(D0_FACSatlas)

#Preprocessing workflow
D0_FACS <- CreateSeuratObject(counts = D0_FACSatlas, project = "MouseAtlas", min.cells = 3, min.features = 200)
D0_FACS
D0_FACS@assays$RNA@data <- D0_FACS@assays$RNA@counts
D0_FACS[["percent.mt"]] <- PercentageFeatureSet(D0_FACS, pattern = "^MT")
head(D0_FACS@meta.data, 5)
VlnPlot(D0_FACS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(D0_FACS, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(D0_FACS, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
D0_FACS <- subset(D0_FACS, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Variable Features
D0_FACS <- FindVariableFeatures(D0_FACS, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(D0_FACS), 10)
plot1 <- VariableFeaturePlot(D0_FACS)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Linear dimension reduction
all.genes <- rownames(D0_FACS)
D0_FACS <- ScaleData(D0_FACS, features = all.genes)
D0_FACS <- RunPCA(D0_FACS, features = VariableFeatures(object = D0_FACS))
print(D0_FACS[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(D0_FACS, dims = 1:2, reduction = "pca")
DimPlot(D0_FACS, reduction = "pca")
DimHeatmap(D0_FACS, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(D0_FACS, dims = 1:15, cells = 500, balanced = TRUE)

#Dimensionality
#D0_FACS <- JackStraw(D0_FACS, num.replicate = 1:10)
#D0_FACS <- ScoreJackStraw(D0_FACS, dims = 1:20)
#JackStrawPlot(D0_FACS, dims = 1:15)
ElbowPlot(D0_FACS)

#Cluster cells
D0_FACS <- FindNeighbors(D0_FACS, dims = 1:10)
D0_FACS <- FindClusters(D0_FACS, resolution = 0.5)
head(Idents(D0_FACS), 5)

#Non-linear dimensional reduction (UMAP/tSNE)
D0_FACS <- RunUMAP(D0_FACS, dims = 1:10)
DimPlot(D0_FACS, reduction = "umap")

#Cluster for cell marking
cluster1.markers <- FindMarkers(D0_FACS, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
cluster5.markers <- FindMarkers(D0_FACS, ident.1 = 5, ident.2 = 0, min.pct = 0.25)
head(cluster5.markers, n = 5)
D0_FACS.markers <- FindAllMarkers(D0_FACS.markers, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
D0_FACS.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
cluster1.markers <- FindMarkers(D0_FACS, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
