library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)

mat_PFC <- read.csv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE124nnn/GSE124952/suppl/GSE124952_expression_matrix.csv.gz")
mat_PFC <- mat_PFC %>%
  as.data.frame() %>%
  column_to_rownames('') %>%
  as.matrix() %>%
  t() 
mat_PFC[1:5, 1:5]
    