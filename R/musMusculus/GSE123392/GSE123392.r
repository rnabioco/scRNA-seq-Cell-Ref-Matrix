library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(digest)
library(here)

proj_dir <- here()
GSE123392 <- file.path(proj_dir, "Reference-Matrix-Generation", "data", "GSE123392_RAW", "Matrices", "GSM3502587_AB2787M.txt.gz")
GSE123392Meta <- file.path(proj_dir, "Reference-Matrix-Generation", "data", "GSE123392_RAW", "Metadata", "GSE123392_scRNAseq_metadata.txt")

mat_GSE123392 <- read.table(GSE123392, sep = " ")
mat_GSE123392 <- as.matrix(mat_GSE123392)
mat_GSE123392[1:5, 1:5]

meta_GSE123392 <- read.table(GSE123392Meta, sep = " ")
meta_GSE123392
ncol(mat_GSE123392)

source("~/Reference-Matrix-Generation/R/utils/utils.r")
checkRawCounts(as.matrix(mat_GSE123392))

GSE123392Normalized <- NormalizeData(mat_GSE123392)

new_ref_matrix <- average_clusters(mat = GSE123392Normalized, metadata = meta_GSE123392$celltype, if_log = FALSE)
head(new_ref_matrix)
tail(new_ref_matrix)
saveRDS(new_ref_matrix, "GSE123392.rds")
