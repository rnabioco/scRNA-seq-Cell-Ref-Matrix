library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)

mat_LiverEndothelial <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129933/suppl/GSE129933_count_matrix.tsv.gz")
mat_LiverEndothelial <- mat_LiverEndothelial %>%
  as.data.frame() %>%
  column_to_rownames('gene') %>%
  as.matrix() %>%
  t()
mat_LiverEndothelial[1:5, 1:5]

meta_LiverEndothelial <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE129nnn/GSE129933/suppl/GSE129933_cell_metadata.tsv.gz")
meta_LiverEndothelial
sum(colnames(mat_LiverEndothelial) %in% meta_LiverEndothelial$cluster)
ncol(mat_LiverEndothelial)

new_ref_matrix <- average_clusters(mat = mat_LiverEndothelial, metadata = meta_LiverEndothelial, if_log = FALSE)
head(new_ref_matrix)
tail(new_ref_matrix)
saveRDS(new_ref_matrix, "GSE129933.rds")
