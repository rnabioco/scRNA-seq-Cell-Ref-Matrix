library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(digest)
library(here)

mat_Lymphocyte <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147668/suppl/GSE147668_table_counts_lungimmune.tsv.gz")
mat_Lymphocyte <- mat_Lymphocyte %>%
  #as.data.frame() %>%
  column_to_rownames('X1')
  #as.matrix() %>%
  #t() 
mat_Lymphocyte[1:5, 1:5]

meta_Lymphocyte <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147668/suppl/GSE147668_table_cellmetadata_lungimmune.tsv.gz")

SeuratLymphocyte <- CreateSeuratObject(counts = mat_Lymphocyte, project = "mat_Lymphocyte", min.cells = 3, min.features = 200)
SeuratLymphocyte

SeuratLymphocyte[["percent.mt"]] <- PercentageFeatureSet(SeuratLymphocyte, pattern = "mt")
# Visualize QC metrics as a violin plot
VlnPlot(SeuratLymphocyte, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(SeuratLymphocyte, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SeuratLymphocyte, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

SeuratLymphocyte <- NormalizeData(SeuratLymphocyte)
SeuratLymphocyte <- FindVariableFeatures(SeuratLymphocyte, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(SeuratLymphocyte), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(SeuratLymphocyte)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(SeuratLymphocyte)
SeuratLymphocyte <- ScaleData(SeuratLymphocyte, features = all.genes)
SeuratLymphocyte <- RunPCA(SeuratLymphocyte, features = VariableFeatures(object = SeuratLymphocyte))
print(SeuratLymphocyte[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(SeuratLymphocyte, dims = 1:2, reduction = "pca")
DimPlot(SeuratLymphocyte, reduction = "pca")
DimHeatmap(SeuratLymphocyte, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(SeuratLymphocyte, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(SeuratLymphocyte)
SeuratLymphocyte <- FindNeighbors(SeuratLymphocyte, dims = 1:10)
SeuratLymphocyte <- FindClusters(SeuratLymphocyte, resolution = 0.5)
head(Idents(SeuratLymphocyte), 5)

SeuratLymphocyte <- RunUMAP(SeuratLymphocyte, dims = 1:10)
DimPlot(SeuratLymphocyte, reduction = "umap")

cluster1.markers <- FindMarkers(SeuratLymphocyte, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
cluster5.markers <- FindMarkers(SeuratLymphocyte, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
SeuratLymphocyte.markers <- FindAllMarkers(SeuratLymphocyte, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SeuratLymphocyte.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
cluster1.markers <- FindMarkers(SeuratLymphocyte, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
top10 <- SeuratLymphocyte.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(SeuratLymphocyte, features = top10$gene) + NoLegend()
saveRDS(SeuratLymphocyte, file = "lymphocyte.rds")

proj_dir <- file.path("Reference-Matrix-Generation", "atlas", "musMusculus", "MouseAtlas.rds")
mouseAtlas <- readRDS(proj_dir)
res <- clustify(
  input = SeuratLymphocyte,          # a Seurat object
  ref_mat = mouseAtlas,         # matrix of RNA-seq expression data for each cell type
  cluster_col = "RNA_snn_res.0.5", # name of column in meta.data containing cell clusters
  obj_out = TRUE               # output Seurat object with cell type inserted as "type" column
)
res@meta.data[1:10, ]
saveRDS(res@meta.data, file = "clustifyLymphocyte.rds")

new.cluster.ids <- res@meta.data$type
names(new.cluster.ids) <- levels(SeuratLymphocyte)
SeuratLung <- RenameIdents(SeuratLymphocyte, new.cluster.ids)
DimPlot(SeuratLymphocyte, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

new.cluster.ids <- meta_Lymphocyte$`Cell Subtype`
names(new.cluster.ids) <- levels(SeuratLymphocyte)
SeuratLung <- RenameIdents(SeuratLymphocyte, new.cluster.ids)
DimPlot(SeuratLymphocyte, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
