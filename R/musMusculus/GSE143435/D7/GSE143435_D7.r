library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)

D7_FACSatlas <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE143nnn/GSE143435/suppl/GSE143435_DeMicheli_D7_FACSatlas_normalizeddata.txt.gz")
D7_FACSatlas <- D7_FACSatlas %>%
  #as.data.frame() %>%
  column_to_rownames('X1')
  #as.matrix() %>%
  #t()
D7_FACSatlas[1:5, 1:5]

D7_FACSatlasMetadata <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE143nnn/GSE143435/suppl/GSE143435_DeMicheli_D7_FACSatlas_metadata.txt.gz")
D7_FACSatlasMetadata
sum(colnames(D7_FACSatlas) %in% D7_FACSatlasMetadata$X1)
ncol(D7_FACSatlas)

source("~/Reference-Matrix-Generation/R/utils/utils.r")
checkRawCounts(as.matrix(D7_FACSatlas))

GSE143435_D7Normalized <- NormalizeData(D7_FACSatlas)
GSE143435_D7Normalized

#Reference matrix build
new_ref_matrix <- average_clusters(mat = GSE143435_D7Normalized, metadata = D7_FACSatlasMetadata$cell_annotation, if_log = TRUE) #Using clustifyr seurat_ref function
head(new_ref_matrix)
tail(new_ref_matrix)            
saveRDS(new_ref_matrix, "GSE143435D7.rds")
