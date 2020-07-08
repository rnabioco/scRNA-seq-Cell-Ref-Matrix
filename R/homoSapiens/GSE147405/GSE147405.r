library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(digest)

KinaseScreen <- read_csv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_OVCA420_TNF_TimeCourse_UMI_matrix.csv.gz")
KinaseScreen <- KinaseScreen %>%  
  #as.data.frame() %>% 
  column_to_rownames('X1')  
  #as.matrix() %>% 
  #t()
KinaseScreen[1:5, 1:5]

meta_KinaseScreen <- read_csv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_OVCA420_EGF_KinaseScreen_metadata.csv.gz")
sum(colnames(KinaseScreen) %in% meta_KinaseScreen$CellLine)
ncol(KinaseScreen)

source("~/Reference-Matrix-Generation/R/utils/check.r")
checkRawCounts(as.matrix(KinaseScreen))

new_ref_matrix <- average_clusters(mat = KinaseScreen, metadata = meta_KinaseScreen$CellLine, if_log = TRUE)
new_ref_matrix_hashed <- average_clusters(mat = KinaseScreen, metadata = meta_KinaseScreen$CellLine, if_log = TRUE)
head(new_ref_matrix)
tail(new_ref_matrix)
newcols <- sapply(colnames(new_ref_matrix_hashed), digest, algo = "sha1")
colnames(new_ref_matrix_hashed) <- newcols
head(new_ref_matrix_hashed)
tail(new_ref_matrix_hashed)
saveRDS(new_ref_matrix_hashed, "GSE147405Hashed.rds")
saveRDS(new_ref_matrix, "GSE147405.rds")
