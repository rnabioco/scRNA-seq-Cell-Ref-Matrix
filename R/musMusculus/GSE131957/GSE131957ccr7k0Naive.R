library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(digest)
library(here)

proj_dir <- here()
GSE131957 <- file.path(proj_dir, "Reference-Matrix-Generation", "data", "GSE131957_Raw 2", "GSM3832739_ccr7ko_naive_gex.csv.gz")

mat_mouseLungLesions <- read.csv(GSE131957)
mat_mouseLungLesions <- mat_mouseLungLesions %>%
  # as.data.frame() %>%
  column_to_rownames('X')
#  as.matrix() %>%
#  t()
mat_mouseLungLesions[1:5, 1:5]

meta_mouseLungLesions <- read_csv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131957/suppl/GSE131957_single_cell_metadata.csv.gz")
meta_mouseLungLesions <- meta_mouseLungLesions %>% column_to_rownames("X1") %>% t() %>% as.data.frame()
meta2 <- meta_mouseLungLesions[colnames(mat_mouseLungLesions),]
meta_mouseLungLesions
sum(colnames(mat_mouseLungLesions) %in% meta_mouseLungLesions$cluster)
ncol(mat_mouseLungLesions)

source("~/Reference-Matrix-Generation/R/utils/utils.r")
checkRawCounts(as.matrix(mat_mouseLungLesions))

GSE131957WTNaiveNormalized <- NormalizeData(mat_mouseLungLesions)

new_ref_matrix <- average_clusters(mat = GSE131957WTNaiveNormalized, metadata = meta2, cluster_col = "cluster_annotation", if_log = FALSE)
new_ref_matrix_hashed <- average_clusters(mat = mat_mouseLungLesions, metadata = meta_mouseLungLesions, if_log = FALSE)
head(new_ref_matrix)
tail(new_ref_matrix)
newcols <- sapply(colnames(new_ref_matrix_hashed), digest, algo = "sha1")
colnames(new_ref_matrix_hashed) <- newcols
head(new_ref_matrix_hashed)
tail(new_ref_matrix_hashed)
saveRDS(new_ref_matrix_hashed, "GSE131957Hashed.rds")
saveRDS(new_ref_matrix, "GSE131957ccr7koNaive.rds")
