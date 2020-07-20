library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(digest)
library(here)

mat_subregions <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113069/suppl/GSE113069_combinedTables.txt.gz")
mat_subregions <- mat_subregions %>%
  # as.data.frame() %>%
  column_to_rownames('B11__d1')
#  as.matrix() %>%
#  t()
mat_subregions[1:5, 1:5]

meta_subregions <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE113nnn/GSE113069/suppl/GSE113069_metadata.txt.gz")
meta_subregions <- meta_subregions %>% mutate(id = str_remove(str_c(id, file, sep = "__"), "\\.dat")) %>% column_to_rownames("id")
meta_subregions2 <- meta_subregions[colnames(mat_subregions),]
meta_subregions2
ncol(mat_subregions)

source("~/Reference-Matrix-Generation/R/utils/utils.r")
checkRawCounts(as.matrix(mat_subregions))

GSE113069Normalized <- NormalizeData(mat_subregions)

new_ref_matrix <- average_clusters(mat = GSE113069Normalized, metadata = meta_subregions2$note2, if_log = FALSE)
head(new_ref_matrix)
tail(new_ref_matrix)
saveRDS(new_ref_matrix, "GSE113069.rds")
