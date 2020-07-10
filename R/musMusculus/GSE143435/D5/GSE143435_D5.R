library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)

D5_FACSatlas <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE143nnn/GSE143435/suppl/GSE143435_DeMicheli_D5_FACSatlas_normalizeddata.txt.gz")
D5_FACSatlas <- D5_FACSatlas %>%
  #as.data.frame() %>%
  column_to_rownames('X1')
  #as.matrix() %>%
  #t()
D5_FACSatlas[1:5, 1:5]

D5_FACSatlasMetadata <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE143nnn/GSE143435/suppl/GSE143435_DeMicheli_D5_FACSatlas_metadata.txt.gz")
D5_FACSatlasMetadata
sum(colnames(D5_FACSatlas) %in% D5_FACSatlasMetadata$X1)
ncol(D5_FACSatlas)

source("~/Reference-Matrix-Generation/R/utils/utils.r")
checkRawCounts(as.matrix(D5_FACSatlas))

GSE143435_D5Normalized <- NormalizeData(D5_FACSatlas)
GSE143435_D5Normalized

#Reference matrix build
new_ref_matrix <- average_clusters(mat = GSE143435_D5Normalized, metadata = D5_FACSatlasMetadata$cell_annotation, if_log = TRUE) #Using clustifyr seurat_ref function
head(new_ref_matrix)
tail(new_ref_matrix)            
saveRDS(new_ref_matrix, "GSE143435D5.rds")
