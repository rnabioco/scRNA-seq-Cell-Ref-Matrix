library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(digest)

mat_Lung <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113049/suppl/GSE113049_count_matrix.tsv.gz", skip = 1, col_names = F)
ids_mat_Lung <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113049/suppl/GSE113049_count_matrix.tsv.gz", n_max = 1, col_names = F)
mat_Lung <- mat_Lung %>%
  #as.data.frame() %>%
  column_to_rownames('X1')
  #as.matrix() %>%
  #t() 
mat_Lung[1:5, 1:5]
colnames(mat_Lung) <- ids_mat_Lung[1,] %>% unlist() %>% as.vector()

meta_Lung <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113049/suppl/GSE113049_cell_metadata.tsv.gz")
meta_Lung
sum(colnames(mat_Lung) %in% meta_Lung$cell_type)
ncol(mat_Lung)

source("~/Reference-Matrix-Generation/R/utils/check.r")
checkRawCounts(as.matrix(mat_Lung))

GSE113049Normalized <- NormalizeData(mat_Lung)
GSE113049Normalized

new_ref_matrix <- average_clusters(mat = GSE113049Normalized, metadata = meta_Lung$cell_type, if_log = F, output_log = F)
new_ref_matrix_hashed <- average_clusters(mat = mat_Lung, metadata = meta_Lung$cell_type, if_log = TRUE)
head(new_ref_matrix)
tail(new_ref_matrix)
newcols <- sapply(colnames(new_ref_matrix_hashed), digest, algo = "sha1")
colnames(new_ref_matrix_hashed) <- newcols
head(new_ref_matrix_hashed)
tail(new_ref_matrix_hashed)
saveRDS(new_ref_matrix_hashed, "GSE113049Hashed.rds")
saveRDS(new_ref_matrix, "GSE113049.rds")
