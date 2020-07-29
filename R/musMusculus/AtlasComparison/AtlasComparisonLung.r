library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(digest)
library(here)
library(Matrix)
library(ggplot2)
library(cowplot)
library(rlang)

proj_dir <- here()
proj_dir <- file.path("Reference-Matrix-Generation", "data", "atlasComparison", "GSE141259")

mat_Lung <- Read10X(data.dir = proj_dir, gene.column = 1)
SeuratLung <- CreateSeuratObject(counts = mat_Lung, project = "mat_Lung", min.cells = 3, min.features = 200)
SeuratLung

meta_Lung <- read_csv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE141nnn/GSE141259/suppl/GSE141259_WholeLung_cellinfo.csv.gz")

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
saveRDS(SeuratLung, file = "lung.rds")

proj_dir <- file.path("Reference-Matrix-Generation", "atlas", "musMusculus", "MouseAtlas.rds")
mouseAtlas <- readRDS(proj_dir)
res <- clustify(
  input = SeuratLung,          # a Seurat object
  ref_mat = mouseAtlas,         # matrix of RNA-seq expression data for each cell type
  cluster_col = "RNA_snn_res.0.5", # name of column in meta.data containing cell clusters
  obj_out = TRUE              # output SCE object with cell type inserted as "type" column
)
res@meta.data[1:10, ]
saveRDS(res@meta.data, file = "clustifyLung.rds")

DimPlot(SeuratLung, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

unique(unlist(res@meta.data$type))

#new.cluster.ids <- c("Macrophage (GSE113049)", "Other Injured AEC2 (GSE113049)", "lung endothelial cell-Lung (Tabula-Muris-drop)", "stromal cell-Lung (Tabula-Muris-drop)", "Ciliated (GSE113049)", "Endothelial/Fibroblast (GSE113049)", "ciliated columnar cell of tracheobronchial tree-Lung (Tabula-Muris-drop)", "T cell-Lung (Tabula-Muris-drop)")
#names(new.cluster.ids) <- levels(SeuratLung)
#SeuratLung <- RenameIdents(SeuratLyung, new.cluster.ids)
Idents(SeuratLung) <- res@meta.data$type
InferredTypes <- DimPlot(SeuratLung, reduction = "umap", label = TRUE, pt.size = 0.5)

unique(unlist(meta_Lung$cell.type))

#new.cluster.ids <- c("NK cell", "T cell", "B cell", "Mac III", "neutrophil", "IL cell", "Mac V", "basophil", "Mac II", "mast cell", "Mac IV", "DC III", "DC II", "DC I", "Mac I")
#names(new.cluster.ids) <- levels(SeuratLung)
#SeuratLung <- RenameIdents(SeuratLung, new.cluster.ids)
Idents(SeuratLung) <- meta_Lung$cell.type
MetadataTypes <- DimPlot(SeuratLung, reduction = "umap", label = TRUE, pt.size = 0.5)

prow <- plot_grid(
  InferredTypes + theme(legend.position="none"),
  MetadataTypes + theme(legend.position="none"),
  align = 'vh',
  labels = c("Lung Inferred", "Lung Metadata"),
  hjust = -1,
  nrow = 1
  )
prow
legendA <- get_legend(
  # create some space to the left of the legend
  InferredTypes + theme(legend.box.margin = margin(0, 0, 0, 12))
)
legendB <- get_legend(
  MetadataTypes + theme(legend.box.margin = margin(0,0,0,12))
)
prow2 <- plot_grid(
  legendA, 
  legendB
  )
plot_grid(prow, prow2, ncol = 1, rel_heights = c(1, .1))
