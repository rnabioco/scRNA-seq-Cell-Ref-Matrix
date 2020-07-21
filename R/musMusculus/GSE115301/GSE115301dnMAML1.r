library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(digest)
library(here)

proj_dir <- here()
GSE115301 <- file.path(proj_dir, "Reference-Matrix-Generation", "data", "GSE115301", "GSE115301_dnMAML1_EV_10x_count.txt")
GSE115301Meta <- file.path(proj_dir, "Reference-Matrix-Generation", "data", "GSE115301", "GSE115301_dnMAML1_EV_10x_metadata.txt")

mat_GSE115301 <- read.table(GSE115301, sep = "\t", row.names = 1)
mat_GSE115301 <- as.matrix(mat_GSE115301)
mat_GSE115301[1:5, 1:5]

meta_GSE115301 <- read.table(GSE115301Meta, sep = "\t", row.names = 1)
meta_GSE115301
ncol(mat_GSE115301)

source("~/Reference-Matrix-Generation/R/utils/utils.r")
checkRawCounts(as.matrix(mat_GSE115301))

GSE115301Normalized <- NormalizeData(mat_GSE115301)

new_ref_matrix <- average_clusters(mat = GSE115301Normalized, metadata = meta_GSE115301$V3, if_log = FALSE)
head(new_ref_matrix)
tail(new_ref_matrix)
saveRDS(new_ref_matrix, "GSE115301.rds")
