library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(digest)
library(here)

proj_dir <- here()
GSE150768 <- file.path(proj_dir, "Reference-Matrix-Generation", "data", "GSE150768", "GSE150768_scRNAseq_rawmatrix.txt")
GSE150768Meta <- file.path(proj_dir, "Reference-Matrix-Generation", "data", "GSE150768", "GSE150768_scRNAseq_metadata.txt")

mat_atheroscleroticlesions <- read.table(GSE150768, sep = " ")
mat_atheroscleroticlesions <- as.matrix(mat_atheroscleroticlesions)
mat_atheroscleroticlesions[1:5, 1:5]

meta_atheroscleroticlesions <- read.table(GSE150768Meta, sep = " ")
meta_atheroscleroticlesions
ncol(mat_atheroscleroticlesions)

source("~/Reference-Matrix-Generation/R/utils/utils.r")
checkRawCounts(as.matrix(mat_atheroscleroticlesions))

GSE150768Normalized <- NormalizeData(mat_atheroscleroticlesions)

new_ref_matrix <- average_clusters(mat = GSE150768Normalized, metadata = meta_atheroscleroticlesions$celltype, if_log = FALSE)
head(new_ref_matrix)
tail(new_ref_matrix)
saveRDS(new_ref_matrix, "GSE150768.rds")
