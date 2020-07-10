library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)

D2_FACSatlas <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE143nnn/GSE143435/suppl/GSE143435_DeMicheli_D2_FACSatlas_normalizeddata.txt.gz")
D2_FACSatlas <- D2_FACSatlas %>%
  #as.data.frame() %>%
  column_to_rownames('X1')
  #as.matrix() %>%
  #t()
D2_FACSatlas[1:5, 1:5]

D2_FACSatlasMetadata <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE143nnn/GSE143435/suppl/GSE143435_DeMicheli_D2_FACSatlas_metadata.txt.gz")
D2_FACSatlasMetadata
sum(colnames(D2_FACSatlas) %in% D2_FACSatlasMetadata$X1)
ncol(D2_FACSatlas)

source("~/Reference-Matrix-Generation/R/utils/utils.r")
checkRawCounts(as.matrix(D2_FACSatlas))

GSE143435_D2Normalized <- NormalizeData(D2_FACSatlas)
GSE143435_D2Normalized

#Reference matrix build
new_ref_matrix <- average_clusters(mat = GSE143435_D2Normalized, metadata = D2_FACSatlasMetadata$cell_annotation, if_log = TRUE) #Using clustifyr seurat_ref function
head(new_ref_matrix)
tail(new_ref_matrix)            
saveRDS(new_ref_matrix, "GSE143435D2.rds")

#Preprocessing workflow
D2_FACS <- CreateSeuratObject(counts = D2_FACSatlas %>% t(), project = "MouseAtlas", min.cells = 3, min.features = 200)
D2_FACS
D2_FACS@assays$RNA@data <- D2_FACS@assays$RNA@counts
D2_FACS[["percent.mt"]] <- D2_FACSatlasMetadata$percent_mito
head(D2_FACS@meta.data, 5)
VlnPlot(D2_FACS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(D2_FACS, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(D2_FACS, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
D2_FACS <- subset(D2_FACS, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Variable Features
D2_FACS <- FindVariableFeatures(D2_FACS, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(D2_FACS), 10)
plot1 <- VariableFeaturePlot(D2_FACS)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Linear dimension reduction
all.genes <- rownames(D2_FACS)
D2_FACS <- ScaleData(D2_FACS, features = all.genes)
D2_FACS <- RunPCA(D2_FACS, features = VariableFeatures(object = D2_FACS))
print(D2_FACS[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(D2_FACS, dims = 1:2, reduction = "pca")
DimPlot(D2_FACS, reduction = "pca")
DimHeatmap(D2_FACS, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(D2_FACS, dims = 1:15, cells = 500, balanced = TRUE)

#Dimensionality
#D2_FACS <- JackStraw(D2_FACS, num.replicate = 1:10)
#D2_FACS <- ScoreJackStraw(D2_FACS, dims = 1:20)
#JackStrawPlot(D2_FACS, dims = 1:15)
ElbowPlot(D2_FACS)

#Cluster cells
D2_FACS <- FindNeighbors(D2_FACS, dims = 1:10)
D2_FACS <- FindClusters(D2_FACS, resolution = 0.5)
head(Idents(D2_FACS), 5)

#Non-linear dimensional reduction (UMAP/tSNE)
D2_FACS <- RunUMAP(D2_FACS, dims = 1:10)
DimPlot(D2_FACS, reduction = "umap")

#Cluster for cell marking
cluster1.markers <- FindMarkers(D2_FACS, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
cluster5.markers <- FindMarkers(D2_FACS, ident.1 = 5, ident.2 = 0, min.pct = 0.25)
head(cluster5.markers, n = 5)
D2_FACS.markers <- FindAllMarkers(D2_FACS.markers, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
D2_FACS.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
cluster1.markers <- FindMarkers(D2_FACS, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

#Assigning cell types to identity to clusters
# get shared cell ids
shared_cell_ids <- intersect(rownames(D2_FACS@meta.data), D2_FACSatlasMetadata$X1)
# subset metadata
D2_FACSatlasMetadata <- filter(D2_FACSatlasMetadata, X1 %in% shared_cell_ids)
# reorder metadata
reorder_idx <- match(rownames(D2_FACS@meta.data), D2_FACSatlasMetadata$X1)
D2_FACSatlasMetadata <- D2_FACSatlasMetadata[reorder_idx, ]
# verify reordering
all(rownames(D2_FACS@meta.data) == D2_FACSatlasMetadata$X1)
# add major_cell_lineage vector to meta.data
D2_FACS@meta.data$annotated <- D2_FACSatlasMetadata$cell_annotation
head(D2_FACS@meta.data$annotated)
new.cluster.ids <- D2_FACS@meta.data$annotated
Idents(D2_FACS) <- "annotated"

#Annotated UMAP
DimPlot(D2_FACS, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
head(Idents(D2_FACS), 5)

#Reference matrix build
new_ref_matrix <- average_clusters(mat = D2_FACSatlas, metadata = D2_FACS@meta.data$annotated, if_log = TRUE) #Using clustifyr seurat_ref function
head(new_ref_matrix)
tail(new_ref_matrix)            
saveRDS(new_ref_matrix, "GSE143435D2.rds")
