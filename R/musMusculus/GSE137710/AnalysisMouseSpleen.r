library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)
library(usethis)

# duplicate cell ids (different samples) in this data, but order in mat and meta are the same
mat_mousespleenDC <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137710/suppl/GSE137710_mouse_spleen_counts_normalized_4464x11755.tsv.gz")
mat_mousespleenDC <- mat_mousespleenDC %>%
  as.data.frame() %>%
  column_to_rownames('X1') %>%
  as.matrix() %>%
  t()
mat_mousespleenDC[1:5,1:5]

meta_mouseSpleenDC <- read_tsv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE137nnn/GSE137710/suppl/GSE137710_mouse_spleen_cell_metadata_4464x9.tsv.gz")
sum(colnames(mat_mousespleenDC) %in% meta_mouseSpleenDC$major_cell_lineage)
ncol(mat_humanSpleenDC)

#Reference matrix build
new_ref_matrix <- average_clusters(mat = mat_mousespleenDC, metadata = meta_mouseSpleenDC$major_cell_lineage, if_log = TRUE) #Using clustifyr seurat_ref function
head(new_ref_matrix)
tail(new_ref_matrix)            
saveRDS(new_ref_matrix, "GSE137710MouseSpleen.rds")