library(dplyr)
library(Seurat)
library(patchwork)
library(clustifyr)
library(tidyverse)

KinaseScreen <- read_csv("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/GSE147405_OVCA420_TNF_TimeCourse_UMI_matrix.csv.gz")
KinaseScreen <- KinaseScreen %>%  
  as.data.frame() %>% 
  column_to_rownames('X1') %>% 
  as.matrix() %>% 
  t()
KinaseScreen[1:5, 1:5]
