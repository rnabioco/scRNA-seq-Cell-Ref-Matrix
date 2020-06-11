library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)

mat_PFC <- read_csv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124952/suppl/GSE124952_expression_matrix.csv.gz")
#Error in read.table(file = file, header = header, sep = sep, quote = quote,  : 
#more columns than column names
#In addition: Warning message:
  #In read.table(file = file, header = header, sep = sep, quote = quote,  :
#                  line 1 appears to contain embedded nulls
mat_PFC <- mat_PFC %>%
  as.data.frame() %>%
  column_to_rownames('X1') %>%
  as.matrix() %>%
  t() 
mat_PFC[1:5, 1:5]

meta_PFC <- read_csv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124952/suppl/GSE124952_meta_data.csv.gz")
meta_PFC
sum(colnames(mat_PFC) %in% meta_PFC$X1)
ncol(mat_PFC)

#Preprocessing workflow
PFC <- CreateSeuratObject(counts = mat_PFC, project = "MousePFC", min.cells = 3, min.features = 200)
PFC
PFC@meta.data
PFC@assays$RNA@data <- PFC@assays$RNA@counts
PFC[["percent.mt"]] <- meta_PFC$percent.mito #Error: Cannot add more or fewer cell meta.data information without values being named with cell names
head(PFC@meta.data, 5)
