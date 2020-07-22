library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(digest)
library(here)

proj_dir <- here()
GSE115746 <- file.path(proj_dir, "Reference-Matrix-Generation", "data", "GSE115746", "GSE115746_cells_exon_counts.csv")
GSE115746Meta <- file.path(proj_dir, "Reference-Matrix-Generation", "data", "GSE115746", "GSE115746_complete_metadata_28706-cells.csv")

mat_GSE115746 <- read_csv(GSE115746) %>% column_to_rownames("X1")
mat_GSE115746 <- as.matrix(mat_GSE115746)
mat_GSE115746[1:5, 1:5]

meta_GSE115746 <- read_csv(GSE115746Meta, skip = 8)
meta_GSE115746 <- as.data.frame(meta_GSE115746)
names(meta_GSE115746) <- colnames(read_csv(GSE115746Meta, n_max = 1))
meta_GSE115746
ncol(mat_GSE115746)

source("~/Reference-Matrix-Generation/R/utils/utils.r")
checkRawCounts(as.matrix(mat_GSE115746))

GSE115746Normalized <- NormalizeData(mat_GSE115746)

new_ref_matrix <- average_clusters(mat = GSE115746Normalized, metadata = meta_GSE115746$cell_cluster, if_log = FALSE)
head(new_ref_matrix)
tail(new_ref_matrix)
saveRDS(new_ref_matrix, "GSE115746.rds")
