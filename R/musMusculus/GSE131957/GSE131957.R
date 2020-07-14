library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(digest)

untar("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE131957&format=file")
mat_mouseLungLesions <- read.csv(UntarMouseLungLesions)
mat_mouseLungLesions <- mat_mouseLungLesions %>%
  # as.data.frame() %>%
  column_to_rownames('gene')
#  as.matrix() %>%
#  t()
mat_mouseLungLesions[1:5, 1:5]

meta_mouseLungLesions <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE131nnn/GSE131057/suppl/GSE131957_single_cell_metadata.csv.gz")
meta_mouseLungLesions
sum(colnames(mat_mouseLungLesions) %in% meta_mouseLungLesions$cluster)
ncol(mat_mouseLungLesions)

source("~/Reference-Matrix-Generation/R/utils/check.r")
checkRawCounts(as.matrix(mat_mouseLungLesions))

new_ref_matrix <- average_clusters(mat = mat_mouseLungLesions, metadata = meta_mouseLungLesions, if_log = FALSE)
new_ref_matrix_hashed <- average_clusters(mat = mat_mouseLungLesions, metadata = meta_mouseLungLesions, if_log = FALSE)
head(new_ref_matrix)
tail(new_ref_matrix)
newcols <- sapply(colnames(new_ref_matrix_hashed), digest, algo = "sha1")
colnames(new_ref_matrix_hashed) <- newcols
head(new_ref_matrix_hashed)
tail(new_ref_matrix_hashed)
saveRDS(new_ref_matrix_hashed, "GSE131957Hashed.rds")
saveRDS(new_ref_matrix, "GSE131957.rds")