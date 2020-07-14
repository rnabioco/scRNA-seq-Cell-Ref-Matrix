library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(digest)

mat_Leukocyte <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146771/suppl/GSE146771_CRC.Leukocyte.10x.TPM.txt.gz")
mat_Leukocyte <- mat_Leukocyte %>%
  # as.data.frame() %>%
  column_to_rownames('gene')
#  as.matrix() %>%
#  t()
mat_Leukocyte[1:5, 1:5]

meta_Leukocyte <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE146nnn/GSE146771/suppl/GSE146771_CRC.Leukocyte.10x.Metadata.txt.gz")
meta_Leukocyte
sum(colnames(mat_Leukocyte) %in% meta_Leukocyte$cluster)
ncol(mat_Leukocyte)

source("~/Reference-Matrix-Generation/R/utils/check.r")
checkRawCounts(as.matrix(mat_Leukocyte))

new_ref_matrix <- average_clusters(mat = mat_Leukocyte, metadata = meta_Leukocyte, if_log = FALSE)
new_ref_matrix_hashed <- average_clusters(mat = mat_Leukocyte, metadata = meta_Leukocyte, if_log = FALSE)
head(new_ref_matrix)
tail(new_ref_matrix)
newcols <- sapply(colnames(new_ref_matrix_hashed), digest, algo = "sha1")
colnames(new_ref_matrix_hashed) <- newcols
head(new_ref_matrix_hashed)
tail(new_ref_matrix_hashed)
saveRDS(new_ref_matrix_hashed, "GSE146771Hashed.rds")
saveRDS(new_ref_matrix, "GSE146771.rds")