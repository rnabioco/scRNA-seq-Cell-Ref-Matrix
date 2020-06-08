library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)

mat_FibroblastNucleiDiabetes <- read_csv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148946/suppl/GSE148946_nuclei_normdata.csv.gz")
mat_FibroblastNucleiDiabetes <- mat_FibroblastNucleiDiabetes %>%
  as.data.frame() %>%
  column_to_rownames('X1') %>%
  as.matrix() %>%
  t()
mat_FibroblastNucleiDiabetes[1:5, 1:5]

meta_FibroblastNucleiDiabetes <- read_csv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE148nnn/GSE148946/suppl/GSE148946_nuclei_metadata.csv.gz")
sum(colnames(mat_FibroblastNucleiDiabetes) %in% meta_FibroblastNucleiDiabetes$celltype)
ncol(mat_FibroblastNucleiDiabetes)

FibroblastNucleiDiabetes <- CreateSeuratObject(counts = mat_FibroblastNucleiDiabetes, project = "FibroblastNucleiDiabetes", min.cells = 3, min.features = 200)
FibroblastNucleiDiabetes
FibroblastNucleiDiabetes@assays$RNA@data <- FibroblastNucleiDiabetes@assays$RNA@counts
FibroblastNucleiDiabetes[["percent.mt"]] <- PercentageFeatureSet(FibroblastNucleiDiabetes, pattern = "^MT-")
head(FibroblastNucleiDiabetes@meta.data, 5)
VlnPlot(FibroblastNucleiDiabetes, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(FibroblastNucleiDiabetes, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(FibroblastNucleiDiabetes, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
FibroblastNucleiDiabetes <- subset(FibroblastNucleiDiabetes, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)