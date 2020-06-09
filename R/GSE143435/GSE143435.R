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
head(D0_FACS@meta.data, 5)
VlnPlot(D0_FACS, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
