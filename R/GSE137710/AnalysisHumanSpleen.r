library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(usethis)

mat_humanSpleenDC <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137710/suppl/GSE137710_human_spleen_counts_normalized_4465x12476.tsv.gz")
mat_humanSpleenDC <- mat_humanSpleenDC %>%
  as.data.frame() %>%
  column_to_rownames('X1') %>%
  as.matrix() %>%
  t()
mat_humanSpleenDC[1:5,1:5]

meta_humanSpleenDC <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137710/suppl/GSE137710_human_spleen_cell_metadata_4465x9.tsv.gz")
sum(colnames(mat_humanSpleenDC) %in% meta_humanSpleenDC$cell_ID)
ncol(mat_humanSpleenDC)

#Preprocessing workflow
humanSpleen <- CreateSeuratObject(counts = mat_humanSpleenDC, project = "HumanSpleen", min.cells = 3, min.features = 200)
humanSpleen
rm(mat_humanSpleenDC)
gc()
humanSpleen@assays$RNA@data <- humanSpleen@assays$RNA@counts
humanSpleen[["percent.mt"]] <- PercentageFeatureSet(humanSpleen, pattern = "^MT")
head(humanSpleen@meta.data, 5)
VlnPlot(humanSpleen, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(humanSpleen, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(humanSpleen, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
humanSpleen <- subset(humanSpleen, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#Variable Features
humanSpleen <- FindVariableFeatures(humanSpleen, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(humanSpleen), 10)
plot1 <- VariableFeaturePlot(humanSpleen)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#Linear Dimension Reduction
all.genes <- rownames(humanSpleen)
humanSpleen <- ScaleData(humanSpleen, features = all.genes)
humanSpleen <- RunPCA(humanSpleen ,features = VariableFeatures(object = humanSpleen))
print(humanSpleen[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(humanSpleen, dims = 1:2, reduction = "pca")
DimPlot(humanSpleen, reduction = "pca")
DimHeatmap(humanSpleen, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(humanSpleen, dims = 1:15, cells = 500, balanced = TRUE)

#Dimensionality
#humanSpleen <- JackStraw(humanSpleen, num.replicate = 100)
#humanSpleen <- ScoreJackStraw(humanSpleen, dims = 1:20)
#JackStrawPlot(humanSpleen, dims = 1:15)
ElbowPlot(humanSpleen)

#Cluster cells
humanSpleen <- FindNeighbors(humanSpleen, dims = 1:10)
humanSpleen <- FindClusters(humanSpleen, resolution = 0.5)
head(Idents(humanSpleen), 5)

#Non-linear dimensional reduction (UMAP/tSNE)
humanSpleen <- RunUMAP(humanSpleen, dims = 1:10)
DimPlot(humanSpleen, reduction = "umap")

cluster1.markers <- FindMarkers(humanSpleen, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
cluster5.markers <- FindMarkers(humanSpleen, ident.1 = 5, ident.2 = 0, min.pct = 0.25)
head(cluster5.markers, n = 5)
humanSpleen.markers <- FindAllMarkers(humanSpleen, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
humanSpleen.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
cluster1.markers <- FindMarkers(humanSpleen, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
VlnPlot(humanSpleen, features = c("IL4I1", "CXCL9"))
VlnPlot(humanSpleen, features = c("IL4I1", "CXCL9"), slot = "counts", log = TRUE)
FeaturePlot(humanSpleen, features = c("IL4I1", "CXCL9", "C15ORF48", "LAMP3", "FSCN1"))
top10 <- humanSpleen.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(humanSpleen, features = top10$gene) + NoLegend()

#Assigning cell types
#get shared cell ids
shared_cell_ids <- intersect(rownames(humanSpleen@meta.data), meta_humanSpleenDC$cell_ID)
#subset meta_humanSpleenDc
meta_humanSpleenDC <- filter(meta_humanSpleenDC, cell_ID %in% shared_cell_ids)
#reorder meta_humanSpleenDC
reorder_idx <- match(rownames(humanSpleen@meta.data), meta_humanSpleenDC$cell_ID)
meta_humanSpleenDC <- meta_humanSpleenDC[reorder_idx, ]
#verify reordering
all(rownames(humanSpleen@meta.data) == meta_humanSpleenDC$cell_ID)
#add major_cell_lineage vector to meta.data
humanSpleen@meta.data$annotated <- meta_humanSpleenDC$cell_type
head(humanSpleen@meta.data$annotated)
new.cluster.ids <- humanSpleen@meta.data$annotated
Idents(humanSpleen) <- "annotated"
DimPlot(humanSpleen, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
head(Idents(humanSpleen), 5)

new_ref_matrix <- average_clusters(mat = mat_humanSpleenDC, metadata = humanSpleen@meta.data$annotated, if_log = TRUE)
head(new_ref_matrix)
tail(new_ref_matrix)
saveRDS(new_ref_matrix, "humanSpleenRefMatrix.rds")
